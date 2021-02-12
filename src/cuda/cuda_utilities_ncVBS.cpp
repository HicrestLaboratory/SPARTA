//CUDA Utilities and system includes
#include <assert.h>
#include <helper_string.h>  // helper for shared functions common to CUDA Samples

// CUDA runtime
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>

// CUDA and CUBLAS functions
#include <helper_functions.h>
#include <helper_cuda.h>

//others
#include<stdio.h>
#include <iostream>
#include <typeinfo>

//includes from SPARTA
#include "cuda_utilities.h"
#include "cuda_utilities_ncVBS.h"
#include "comp_mats.h"
#include "sparse_utilities.h"

/*ONLY WORKS WITH 
    DataT = float, double, int;
    intT = int
    unsigned_intT = unsigned int;
*/


void cublas_ncVBS_multiply(ncVBS& vbmatA, const DataT* B, int B_cols, int B_lead_dim, DataT* C, int C_lead_dim, float* times, int n_streams_mult, int n_streams_cpy)
{
    cudaDataType_t cuda_type;
    if (typeid(DataT) == typeid(float))    cuda_type = CUDA_R_32F;
    else if (typeid(DataT) == typeid(double))   cuda_type = CUDA_R_64F;
    else if (typeid(DataT) == typeid(int))  cuda_type = CUDA_R_8I;

    cublasGemmAlgo_t cuda_algo = CUBLAS_GEMM_DEFAULT_TENSOR_OP;

    intT A_rows = vbmatA.rows;
    intT A_cols = vbmatA.cols();

    intT B_rows = A_cols;

    intT C_rows = A_rows;
    intT C_cols = B_cols;

    const DataT alpha = 1.0f;
    const DataT beta = 1.0f;

    cublasHandle_t handle;
    checkCudaErrors(cublasCreate(&handle));


    int number_of_time_measures = 6;
    int t_counter = 0;
    //define events
    cudaEvent_t t_measures[number_of_time_measures];

    for (int i = 0; i < number_of_time_measures; i++)
    {
        cudaEventCreate(t_measures + i);
    }

    cudaEventRecord(t_measures[0], 0);
    t_counter++;

    //starting the streams
    n_streams_mult = std::min(n_streams_mult, vbmatA.block_cols);
    cudaStream_t streams_mult[n_streams_mult];
    cudaStream_t streams_cpy[n_streams_cpy];

    for (intT ib = 0; ib < n_streams_mult; ib++)    cudaStreamCreate(&(streams_mult[ib]));

    for (intT ib = 0; ib < n_streams_cpy; ib++)     cudaStreamCreate(&(streams_cpy[ib]));


    //event 1: streams created

    cudaEventRecord(t_measures[t_counter], 0);
    cudaEventSynchronize(t_measures[t_counter]);
    cudaEventElapsedTime(times + t_counter, t_measures[0], t_measures[t_counter]);
    cudaEventDestroy(t_measures[t_counter]);
    t_counter++;


    //allocate d_C, initialize to 0;
    unsigned int size_C_whole = C_rows * C_cols;
    unsigned int mem_size_C_whole = sizeof(DataT) * size_C_whole;
    DataT* d_C_whole;
    checkCudaErrors(cudaMalloc((void**)&d_C_whole, mem_size_C_whole));
    cudaMemset(d_C_whole, 0, mem_size_C_whole);
    
    
    DataT** d_A_blocks = new DataT * [vbmatA.block_cols];
    DataT** d_B_blocks = new DataT * [vbmatA.block_cols];


    //allocate and fill d_B. Set pointers in d_B_blocks to right pos in d_B
    unsigned int size_B = B_rows * B_cols;
    unsigned int mem_size_B = sizeof(DataT) * size_B;
    DataT* d_B;
    checkCudaErrors(cudaMalloc((void**)&d_B, mem_size_B));
    checkCudaErrors(cublasSetVector(
        B_cols * B_rows, sizeof(DataT), B, 1, d_B, 1));
    for (int jb = 0; jb < vbmatA.block_cols; jb++)
    {
        d_B_blocks[jb] = d_B + B_cols* vbmatA.block_width(jb);
    }

    
    //allocate and fill d_A. Set pointers in d_A_blocks to right pos in d_A
    unsigned int size_A_total = vbmatA.tot_elems();
    unsigned int mem_size_A = sizeof(DataT) * size_A_total;
    DataT* d_A;
    checkCudaErrors(cudaMalloc((void**)&d_A, mem_size_A));
    
    DataT* pointerA = d_A;
    intT block_size;
    for (int jb = 0; jb < vbmatA.block_cols; jb++)
    {
        d_A_blocks[jb] = pointerA;
        block_size = (vbmat.col_part[jb + 1] - vbmat.col_part[jb]) * vbmat.nzcount[jb];

        checkCudaErrors(cublasSetVector(
            block_size, sizeof(DataT), vbmatA.mab[jb], 1, d_A_blocks[jb], 1));

        pointerA += block_size;
    }


    //allocate d_C_blocks
    DataT** d_C_blocks = new DataT * [vbmatA.block_cols];
    unsigned int C_block_size;
    unsigned int C_block_mem;
    for (int jb = 0; ib < vbmatA.block_cols; jb++)
    {
        C_block_size = vbmatA.nzcount[jb] * B_cols;
        C_block_mem = sizeof(DataT) * C_block_size;
        checkCudaErrors(cudaMalloc((void**)&d_C_blocks[jb], C_block_mem));
    }



    //event 2: data transfer completed
    cudaEventRecord(t_measures[t_counter], 0);
    cudaEventSynchronize(t_measures[t_counter]);
    cudaEventElapsedTime(times + t_counter, t_measures[0], t_measures[t_counter]);
    cudaEventDestroy(t_measures[t_counter]);
    t_counter++;



    
    for (intT jb = 0; jb < vbmatA.block_cols; jb++)
    {

        int rows_number = vbmatA.nzcount[jb];
        int* rows_indices = vbmatA.nzindex[jb];
        int column_start = vbmatA.col_part[jb];
        int column_end = vbmatA.col_part[jb + 1];
        int column_block_size = column_end - column_start;


        int stream_n = jb % n_streams_mult;
        cublasSetStream(handle, streams_mult[stream_n]);

        checkCudaErrors(
            cublasGemmEx(
                handle, CUBLAS_OP_N, CUBLAS_OP_N,
                B_cols, rows_number, column_block_size,           //m, n, k <-- block_B: m*k   block_A: k*n   block_C: n*m
                &alpha,
                d_B_blocks[jb],                                      // blockB device pointer,
                cuda_type,                                      // blockB datatype
                B_cols,                                  // blockB leading dimension
                d_A_blocks[jb],                                      // blockA device pointer
                cuda_type,                                      // blockA datatype
                column_block_size,                                         // leading dimension
                &beta,
                d_C_blocks[jb], cuda_type,                           // blockC device pointer, blockC type
                B_cols,
                cuda_type,                                      // compute_type
                cuda_algo)
        ); 
      
    }

    cudaDeviceSynchronize();
    //event 3: multiplication end
    cudaEventRecord(t_measures[t_counter], 0);
    cudaEventSynchronize(t_measures[t_counter]);
    cudaEventElapsedTime(times + t_counter, t_measures[0], t_measures[t_counter]);
    cudaEventDestroy(t_measures[t_counter]);
    t_counter++;


    //accumulate onto C_whole, row_by_row
    for (int jb = 0; jb < vbmatA.block_cols; jb++)
    {
        unsigned int rows_number = vbmatA.nzcount[jb];

        for (int i = 0; i < rows_number; i++)
        {
            int stream_n = i % n_streams_cpy;
            cublasSetStream(handle, streams_cpy[stream_n]);

            int C_pos = vbmatA.nzindex[jb][i] * C_cols;
            //TODO make it work for other datatypes
            checkCudaErrors(
                cublasSaxpy(handle, C_cols,
                    &alpha,
                    d_C_blocks[jb] + i * C_cols, 1,
                    d_C_whole + C_pos, 1)
            );
        }
    }
    cublasSetStream(handle, 0);
    cudaDeviceSynchronize();

    //event 4: accumulation end
    cudaEventRecord(t_measures[t_counter], 0);
    cudaEventSynchronize(t_measures[t_counter]);
    cudaEventElapsedTime(times + t_counter, t_measures[0], t_measures[t_counter]);
    cudaEventDestroy(t_measures[t_counter]);
    t_counter++;

    //copy into C
    checkCudaErrors(
        cublasGetVector(
            C_rows*C_cols, sizeof(DataT),
            d_C_whole, 1,
            C, 1
        )
    );


    //event 5: final copy end
    cudaEventRecord(t_measures[t_counter], 0);
    cudaEventSynchronize(t_measures[t_counter]);
    cudaEventElapsedTime(times + t_counter, t_measures[0], t_measures[t_counter]);
    cudaEventDestroy(t_measures[t_counter]);
    t_counter++;



    //cleanup
    cudaEventDestroy(t_measures[0]);

    for (int jb = 0; jb < vbmatA.block_cols; jb++)
    {
        checkCudaErrors(cudaFree(d_C[jb]));
        checkCudaErrors(cudaFree(d_A[jb]));
        checkCudaErrors(cudaFree(d_B[jb]));
    }
    checkCudaErrors(cudaFree(d_C_whole));

    checkCudaErrors(cublasDestroy(handle));

}
