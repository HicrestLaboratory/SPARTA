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


void cublas_ncVBS_multiply(const ncVBS& vbmatA, DataT* B, int B_cols, int B_lead_dim, DataT* C, int C_lead_dim, float& dt, int n_streams_mult = 32, int n_streams_cpy = 32)
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


    n_stream_mult = std::min(n_streams_mult, vbmatA.block_cols);

    for (intT ib = 0; ib < n_streams_mult; ib++)
    {
        cudaStreamCreate(&(streams[ib]));
    }

    for (intT jb = 0; jb < vbmat.block_cols; jb++)
    {
        intT rows_number = vbmat.nzcount[jb];
        intT* rows_indices = vbmat.nzindex[jb];
        intT column_start = vbmat.col_part[jb];
        intT column_end = vbmat.col_part[jb + 1];
        intT column_block_size = column_end - column_start;


        intT stream_n = jb % n_streams_mult
            ;
        cublasSetStream(handle, streams[ib]);

    }



}