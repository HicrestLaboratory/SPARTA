// Utilities and system includes
#include <assert.h>
#include <helper_string.h>  // helper for shared functions common to CUDA Samples

// CUDA runtime
#include <cuda_runtime.h>
#include <cublas_v2.h>

// CUDA and CUBLAS functions
#include <helper_functions.h>
#include <helper_cuda.h>

#include<stdio.h>


#include "cuda_utilities.h"
#include "globheads.h"
#include "protos.h"
#include "utilities.h"

#include <iostream>


using namespace std;


void cublas_blockmat_multiply(const VBSparMat &VBMat, float *X, int X_cols, float *Y){
    int N = VBMat.n, *bsz = VBMat.bsz;
    int block_cols,block_rows,col;
    int mat_n = bsz[N];
    int X_rows = mat_n;
    int Y_rows = mat_n;
    int Y_cols = X_cols;

    const float alpha = 1.0f;
    const float beta = 1.0f;
   
    //loop vertically through block rows
    for(int i = 0; i < N; i++ ) {
        block_rows = bsz[i+1] - bsz[i];


	//allocate device memory for block d_Y
	//-----------------------------------------------
	unsigned int size_dY = block_rows * Y_cols;
	unsigned int mem_size_dY = sizeof(float) * size_dY;

	float *d_Y;
	checkCudaErrors(cudaMalloc((void **) &d_Y, mem_size_dY));

        //-----------------------------------------------


	//loop horizontaly through block columns
        for(int j = 0; j<VBMat.nzcount[i]; j++){
		col = VBMat.ja[i][j];
		block_cols = bsz[col+1] - bsz[col];
		//multiply the block matrices
		 //define the sub-matrices
		const float* block = (VBMat.ba)[i][j];  //access block i,j (stored in column major order).
		const float* blockX = X + bsz[col];     //access the block of X that is going to be multiplied with the (i,j)block of VBMat

		cublas_gemm_custom(block, block_cols, block_rows, block_rows,
                   blockX, X_cols, X_rows,
                   d_Y, block_rows,
		   alpha, beta);
	}

	//retrieve matrix from device and free memory
	/*--------------------------------------*/

	float* blockY = Y + bsz[i];             //i indicates the vertical block of Y that is going to be updated               
	checkCudaErrors(cublasGetMatrix(
                                    block_rows, Y_cols, sizeof(float), d_Y, block_rows, blockY, Y_rows));
   
	checkCudaErrors(cudaFree(d_Y));
	/*-------------------------------------*/

    }
}



//Matrix-Matrix multiplication with cublas. A,B,C are in column-major order.
//Matrix A and B are in host
//Matrix d_C is in device to allow for accumulation of results
int cublas_gemm_custom(const float *A, unsigned int A_cols, unsigned int A_rows, unsigned int lda,
	const float *B, unsigned int B_cols, unsigned int ldb,
	float *d_C, unsigned int ldc,
	const float alpha,
	const float beta)
{
    int block_size = 16;
    cublasStatus_t stat;
    
    //deduce matrices dimensions
    unsigned int B_rows = A_cols;
    unsigned int C_rows = A_rows;
    unsigned int C_cols = B_cols;
    
    unsigned int size_A = A_rows * A_cols;
    unsigned int mem_size_A = sizeof(float) * size_A;

    unsigned int size_B = B_rows * B_cols;
    unsigned int mem_size_B = sizeof(float) * size_B;    
    
    // allocate device memory
    float *d_A, *d_B;
    checkCudaErrors(cudaMalloc((void **) &d_A, mem_size_A));
    checkCudaErrors(cudaMalloc((void **) &d_B, mem_size_B)); 
    
    //copy matrices to device
    checkCudaErrors(cublasSetMatrix(
                                    A_rows, A_cols, sizeof(float), A, lda, d_A, A_rows));
    checkCudaErrors(cublasSetMatrix(
                                    B_rows, B_cols, sizeof(float), B, ldb, d_B, B_rows));

    // setup execution parameters
    dim3 threads(block_size, block_size);
    dim3 grid(C_cols / threads.x, C_rows / threads.y);

    // CUBLAS version 2.0
    cublasHandle_t handle;

    checkCudaErrors(cublasCreate(&handle));

    //Perform gemm operation with cublas
    checkCudaErrors(cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                                A_rows, B_cols, A_cols,
                                &alpha,
                                d_A, A_rows,
                                d_B, B_rows,
                                &beta,
                                d_C, C_rows));
    // copy result from device to host 

    // Destroy the handle
    checkCudaErrors(cublasDestroy(handle));

    // clean up memory
    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_B));
	
    return 0;
}


void randomInit(float *data, int size)
{
    for (int i = 0; i < size; ++i)
        data[i] = rand() / (float)RAND_MAX;
}

