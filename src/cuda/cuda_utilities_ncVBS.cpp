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


void cublas_ncVBS_multiply(ncVBS& vbmatA, DataT* B, int B_cols, int B_lead_dim, DataT* C, int C_lead_dim, float& dt, int n_streams_mult = 32, int n_streams_cpy = 32)
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

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);


    n_streams_mult = std::min(n_streams_mult, vbmatA.block_cols);

    cudaStream_t streams_mult[n_streams_mult];
    cudaStream_t streams_cpy[n_streams_cpy];
    
    for (intT ib = 0; ib < n_streams_mult; ib++)
    {
        cudaStreamCreate(&(streams_mult[ib]));
    }

    for (intT ib = 0; ib < n_streams_cpy; ib++)
    {
        cudaStreamCreate(&(streams_cpy[ib]));
    }

    unsigned int size_C_whole = C_rows * C_cols;
    unsigned int mem_size_C_whole = sizeof(DataT) * size_C_whole;
    DataT* d_C_whole;
    checkCudaErrors(cudaMalloc((void**)&d_C_whole, mem_size_C_whole));
    cudaMemset(d_C_whole, 0, mem_size_C_whole);
    
    
    DataT** d_A = new DataT*[vbmatA.block_cols];
    DataT** d_B = new DataT*[vbmatA.block_cols];
    DataT** d_C = new DataT*[vbmatA.block_cols];


    for (intT jb = 0; jb < vbmatA.block_cols; jb++)
    {

        int rows_number = vbmatA.nzcount[jb];
        int* rows_indices = vbmatA.nzindex[jb];
        int column_start = vbmatA.col_part[jb];
        int column_end = vbmatA.col_part[jb + 1];
        int column_block_size = column_end - column_start;


        int stream_n = jb % n_streams_mult;
        cublasSetStream(handle, streams_mult[stream_n]);


        unsigned int size_A_block = column_block_size * rows_number;
        unsigned int mem_size_A = sizeof(DataT) * size_A_block;

        unsigned int size_B_block = column_block_size * B_cols;
        unsigned int mem_size_B = sizeof(DataT) * size_B_block;

        unsigned int size_C_block = rows_number*B_cols;
        unsigned int mem_size_C = sizeof(DataT) * size_C_block;

        checkCudaErrors(cudaMalloc((void**)&d_A[jb], mem_size_A));
        checkCudaErrors(cudaMalloc((void**)&d_B[jb], mem_size_B));
        checkCudaErrors(cudaMalloc((void**)&d_C[jb], mem_size_C));

        checkCudaErrors(cublasSetVector(
            size_A_block, sizeof(DataT), vbmatA.mab[jb], 1, d_A[jb], 1));

        int B_start_position = column_start * B_cols;
        checkCudaErrors(cublasSetVector(
            size_B_block, sizeof(DataT), B + B_start_position, 1, d_B[jb], 1));

        checkCudaErrors(
            cublasGemmEx(
                handle, CUBLAS_OP_N, CUBLAS_OP_N,
                B_cols, rows_number, column_block_size,           //m, n, k <-- block_B: m*k   block_A: k*n   block_C: n*m
                &alpha,
                d_B[jb],                                      // blockB device pointer,
                cuda_type,                                      // blockB datatype
                B_cols,                                  // blockB leading dimension
                d_A[jb],                                      // blockA device pointer
                cuda_type,                                      // blockA datatype
                column_block_size,                                         // leading dimension
                &beta,
                d_C[jb], cuda_type,                           // blockC device pointer, blockC type
                B_cols,
                cuda_type,                                      // compute_type
                cuda_algo)
        ); 
      
    }
    cudaDeviceSynchronize();

    //accumulate onto C_whole, row_by_row
    for (int jb = 0; jb < vbmatA.block_cols; jb++)
    {
        unsigned int rows_number = vbmatA.nzcount[jb];

        for (int i = 0; i < rows_number; i++)
        {
            int stream_n = i % n_streams_cpy;
            cublasSetStream(handle, streams_cpy[stream_n]);

            unsigned int C_pos = vbmatA.nzindex[jb][i] * C_cols;
            //TODO make it work for other datatypes
            checkCudaErrors(
                cublasSaxpy(handle, C_cols,
                    1.,
                    d_C[jb] + i * C_cols, 1,
                    d_C_whole + C_pos, 1)
            );
        }
    }

    //copy into C
    checkCudaErrors(
        cublasGetMatrix(
            C_rows, C_cols, sizeof(DataT),
            d_C_whole, C_cols, 
            C, C_lead_dim
        )
    );


    cudaDeviceSynchronize();

    cudaEventRecord(stop, 0);

    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&dt, start, stop);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    checkCudaErrors(cudaFree(d_C));
    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_B));

    checkCudaErrors(cublasDestroy(handle));



}