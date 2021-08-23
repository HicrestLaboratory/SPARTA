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
#include "comp_mats.h"
#include "sparse_utilities.h"

/*ONLY WORKS WITH 
    DataT = float, double, int;
    intT = int
*/

void cublas_blockmat_multiply(const VBS &vbmatA, DataT *B, int B_cols, int B_lead_dim, DataT *C, int C_lead_dim, float &dt, int n_streams = 16){
    //multiplies a VBS matrix (vbmatA) and a dense matrix (B); stores into (C)
    //vbmatA:       column-major entries storage;
    //              column-major block_storage; 
    //B:            column-major storage; TODO: allow general storage format (implement through cublas transpose)
    //C:            column-major storage; TODO: allow general storage format (implement through cublas transpose)


    cudaDataType_t cuda_type;
    if (typeid(DataT) == typeid(float))    cuda_type = CUDA_R_32F;
    else if (typeid(DataT) == typeid(double))   cuda_type = CUDA_R_64F;
    else if (typeid(DataT) == typeid(int))  cuda_type = CUDA_R_8I;

    cublasGemmAlgo_t cuda_algo = CUBLAS_GEMM_DEFAULT_TENSOR_OP;

    int A_rows = vbmatA.row_part[vbmatA.block_rows];
    int A_cols = vbmatA.col_part[vbmatA.block_cols];

    intT mat_idx = 0; //keeps writing position for mat
    intT vbmat_idx = 0; //keeps reading position for vbmat 
    intT ja_count = 0; //keeps total nonzero blocks count;

    int B_rows = A_cols;

    int C_rows = A_rows;
    int C_cols = B_cols;

    const DataT alpha = 1.0f;
    const DataT beta = 1.0f;
   
    intT tot_nonzero_blocks = 0; //index for total nonzero blocks

    int rows_in_block, cols_in_block;
    int size_block, mem_size_block;

    //TODO: allocate memory on device
    intT size_A = vbmatA.nztot; //total nonzero entries in vbmat
    intT mem_size_A = sizeof(DataT) * size_A;

    intT size_B = B_rows * B_cols;
    intT mem_size_B = sizeof(DataT) * size_B;

    intT size_C = C_rows * C_cols;
    intT mem_size_C = sizeof(DataT) * size_C;

    cublasHandle_t handle;

    checkCudaErrors(cublasCreate(&handle));

    DataT* d_A, * d_B, * d_C;
    checkCudaErrors(cudaMalloc((void**)&d_A, mem_size_A));
    checkCudaErrors(cudaMalloc((void**)&d_B, mem_size_B));
    checkCudaErrors(cudaMalloc((void**)&d_C, mem_size_C));

    //copy to device the vbmat matrix (nonzero blocks are stored consecutively and in column major format)
    checkCudaErrors(cublasSetVector(
        size_A, sizeof(DataT), vbmatA.mab, 1, d_A, 1));

    //copy B to device (maybe more efficient to copy it block by block?)
    checkCudaErrors(cublasSetMatrix(
        B_rows, B_cols, sizeof(DataT), B, B_lead_dim, d_B, B_rows));

    //initialize cuda events
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);

    //creates streams. Each block rows is assigned a different stream.
    

    if (n_streams > vbmatA.block_rows) n_streams = vbmatA.block_rows;
    cudaStream_t streams[n_streams];
    for (intT ib = 0; ib < n_streams; ib++)
    {
        cudaStreamCreate(&(streams[ib]));
    }

    //loop through all blocks
    for(intT jb = 0; jb < vbmatA.block_cols; jb++ )      //loop horizontally through block columns
    {
        cols_in_block = vbmatA.col_part[jb+1] - vbmatA.col_part[jb];
        const DataT* d_B_block = d_B + vbmatA.col_part[jb];    //access the block of B that is going to be multiplied with blocks of A in column jb

        for(intT nzs = 0; nzs < vbmatA.nzcount[jb]; nzs++)        //loop vertically through nonzero blocks

        {

            intT ib = vbmatA.jab[tot_nonzero_blocks];             //the block row position of a nonzero block 
            tot_nonzero_blocks += 1;
            rows_in_block = vbmatA.row_part[ib + 1] - vbmatA.row_part[ib]; //the row height of the block

            cublasSetStream(handle, streams[ib%n_streams]);               //each block row works on a different stream

            //define the sub-matrices
	        const DataT* d_A_block = d_A + vbmat_idx;           //access the block on d_A.
            DataT* d_C_block = d_C + vbmatA.row_part[ib] ;      //access the block on d_C.




            //multiply the blocks, store result in d_C_block
            checkCudaErrors(
                cublasGemmEx(
                    handle, CUBLAS_OP_N, CUBLAS_OP_N,
                    rows_in_block, B_cols, cols_in_block,           //m, n, k <-- block_A: m*k   block_B: k*n   block_C: m*n
                    &alpha,
                    d_A_block,                                      // blockA device pointer,
                    cuda_type,                                      // blockA datatype
                    rows_in_block,                                  // blockA leading dimension
                    d_B_block,                                      // blockB device pointer
                    cuda_type,                                      // blockB datatype
                    B_rows,                                         // leading dimension
                    &beta,
                    d_C_block, cuda_type,                           // blockC device pointer, blockC type
                    C_rows,
                    cuda_type,                                      // compute_type
                    cuda_algo)
            );                                       
            vbmat_idx += rows_in_block * cols_in_block;
	    }

    }

    //record the elapsed time onto dt
    cudaDeviceSynchronize();
    cudaEventRecord(stop, 0);

    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&dt, start, stop);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);



    //let each stream copy the relevant C block from device
    int stream_id;
    for (int ib = 0; ib < vbmatA.block_rows; ib++)
    {
        stream_id = ib % n_streams;
        cublasSetStream(handle, streams[stream_id]);
        rows_in_block = vbmatA.row_part[ib + 1] - vbmatA.row_part[ib];
        checkCudaErrors(cublasGetMatrix(
            rows_in_block, C_cols, sizeof(DataT), d_C + vbmatA.row_part[ib], C_rows, C + vbmatA.row_part[ib], C_lead_dim));

    }

    cudaDeviceSynchronize();

    checkCudaErrors(cudaFree(d_C));
    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_B));

    checkCudaErrors(cublasDestroy(handle));

}

//Matrix-Matrix multiplication with cublas. A,B,C are in column-major order.
//Matrix A and B are in host
//Matrix d_C is in device to allow for accumulation of results
int cublas_gemm_custom(const DataT *A, unsigned int A_rows, unsigned int A_cols, unsigned int lda,
	const DataT *B, unsigned int B_cols, unsigned int ldb,
	DataT *C, unsigned int ldc,
	const DataT alpha,
	const DataT beta,
    float& dt)
{

    cudaDataType_t cuda_type;
    if (typeid(DataT) == typeid(float))    cuda_type = CUDA_R_32F;
    else if (typeid(DataT) == typeid(double))   cuda_type = CUDA_R_64F;
    else if (typeid(DataT) == typeid(int))  cuda_type = CUDA_R_8I;

    cublasGemmAlgo_t cuda_algo = CUBLAS_GEMM_DEFAULT_TENSOR_OP;

    int mem = 0;
    
    //deduce matrices dimensions
    unsigned int B_rows = A_cols;
    unsigned int C_rows = A_rows;
    unsigned int C_cols = B_cols;

    //allocate memory on device
    //-------------------------------------------------------
    unsigned int size_A = A_rows * A_cols;
    unsigned int mem_size_A = sizeof(DataT) * size_A;

    unsigned int size_B = B_rows * B_cols;
    unsigned int mem_size_B = sizeof(DataT) * size_B;    
  
    unsigned int size_C = C_rows * C_cols;
    unsigned int mem_size_C = sizeof(DataT) * size_C;




    DataT *d_A, *d_B, *d_C;
    checkCudaErrors(cudaMalloc((void **) &d_A, mem_size_A));
    


    checkCudaErrors(cudaMalloc((void **) &d_B, mem_size_B)); 




    checkCudaErrors(cudaMalloc((void **) &d_C, mem_size_C));
    //-------------------------------------------------------



    //copy matrices to device
    checkCudaErrors(cublasSetMatrix(
                                    A_rows, A_cols, sizeof(DataT), A, lda, d_A, A_rows));
    checkCudaErrors(cublasSetMatrix(
                                    B_rows, B_cols, sizeof(DataT), B, ldb, d_B, B_rows));



    // CUBLAS version 2.0
    cublasHandle_t handle;

    checkCudaErrors(cublasCreate(&handle));




    //initialize cuda events
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);



    //Perform gemm operation with cublas
    checkCudaErrors(
        cublasGemmEx(
            handle, CUBLAS_OP_N, CUBLAS_OP_N,
            A_rows, B_cols, A_cols,
            &alpha,
            d_A, cuda_type, A_rows,
            d_B, cuda_type, B_rows,
            &beta,
            d_C, cuda_type, C_rows,
            cuda_type,
            cuda_algo
            ));

    //record the elapsed time onto dt
    cudaDeviceSynchronize();
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&dt, start, stop);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);



    // copy result from device to host 
    checkCudaErrors(cublasGetMatrix(C_rows, C_cols, sizeof(DataT), d_C, C_rows, C, C_rows));



    // clean up memory
    checkCudaErrors(cudaFree(d_C));
    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_B));



    // Destroy the handle
    checkCudaErrors(cublasDestroy(handle));
	
    return 0;
}


int cusparse_gemm_custom(int rows, int cols, int nnz, int* csrRowPtr, int* csrColInd, float* csrVal, float* B, int B_cols, int B_lead_dim, float* C, int C_lead_dim, 
    const float alpha,
    const float beta,
    float& dt)
{


    //allocate memory on device
    unsigned int mem_size_csrVal = sizeof(float) * nnz;
    unsigned int mem_size_csrColInd = sizeof(int) * nnz;
    unsigned int mem_size_csrRowPtr = sizeof(int) * (rows + 1);

    unsigned int B_rows = cols;
    unsigned int size_B = B_rows * B_cols;
    unsigned int mem_size_B = sizeof(float) * size_B;
    
    unsigned int C_rows = rows;
    unsigned int C_cols = B_cols;
    unsigned int size_C = C_rows * C_cols;
    unsigned int mem_size_C = sizeof(float) * size_C;

    // allocate device memory
    int* d_RowPtr, * d_ColInd;
    float* d_Val;

    checkCudaErrors(cudaMalloc((void**)&d_RowPtr, mem_size_csrRowPtr));
    checkCudaErrors(cudaMalloc((void**)&d_ColInd, mem_size_csrColInd));
    checkCudaErrors(cudaMalloc((void**)&d_Val, mem_size_csrVal));

    float * d_B, * d_C;
    checkCudaErrors(cudaMalloc((void**)&d_B, mem_size_B));
    checkCudaErrors(cudaMalloc((void**)&d_C, mem_size_C));

    //copy arrays and matrices to device
    checkCudaErrors(cublasSetVector(
        nnz, sizeof(float),
        csrVal, 1, d_Val, 1));

    checkCudaErrors(cublasSetVector(
        nnz, sizeof(int),
        csrColInd, 1, d_ColInd, 1));

    checkCudaErrors(cublasSetVector(
        (rows + 1), sizeof(int),
        csrRowPtr, 1, d_RowPtr, 1));

    checkCudaErrors(cublasSetMatrix(
        B_rows, B_cols, sizeof(float), B, B_lead_dim, d_B, B_rows));


    if (beta != 0)
    {
        checkCudaErrors(cublasSetMatrix(
            C_rows, C_cols, sizeof(float), C, C_lead_dim, d_C, C_rows));
    }

    cusparseHandle_t handle;
    cusparseSpMatDescr_t matA;
    cusparseDnMatDescr_t matB, matC;

    checkCudaErrors(cusparseCreate(&handle));

    checkCudaErrors(
        cusparseCreateCsr(
            &matA,
            rows,
            cols,
            nnz,
            d_RowPtr,
            d_ColInd,
            d_Val,
            CUSPARSE_INDEX_32I,
            CUSPARSE_INDEX_32I,
            CUSPARSE_INDEX_BASE_ZERO,
            CUDA_R_32F
            )
        );

    checkCudaErrors(
        cusparseCreateDnMat(&matB, cols, B_cols, B_lead_dim, d_B,
            CUDA_R_32F, CUSPARSE_ORDER_COL));

    checkCudaErrors(
        cusparseCreateDnMat(&matC, rows, B_cols, C_lead_dim, d_C,
            CUDA_R_32F, CUSPARSE_ORDER_COL));

    //initialize cuda events
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    
    
    checkCudaErrors(cusparseSpMM(handle,
        CUSPARSE_OPERATION_NON_TRANSPOSE,
        CUSPARSE_OPERATION_NON_TRANSPOSE,
        &alpha, matA, matB, &beta, matC, CUDA_R_32F,
        CUSPARSE_SPMM_ALG_DEFAULT, dBuffer));

    //record the elapsed time onto dt
    cudaDeviceSynchronize();
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&dt, start, stop);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    // destroy matrix/vector descriptors
    CHECK_CUSPARSE(cusparseDestroySpMat(matA));
    CHECK_CUSPARSE(cusparseDestroyDnMat(matB));
    CHECK_CUSPARSE(cusparseDestroyDnMat(matC));
    CHECK_CUSPARSE(cusparseDestroy(handle));

    // copy result from device to host 
    checkCudaErrors(cublasGetMatrix(C_rows, C_cols, sizeof(float), d_C, C_rows, C, C_rows));

    // clean up memor;y
    checkCudaErrors(cudaFree(d_B));
    checkCudaErrors(cudaFree(d_C));
    checkCudaErrors(cudaFree(d_Val));
    checkCudaErrors(cudaFree(d_RowPtr));
    checkCudaErrors(cudaFree(d_ColInd));

    // Destroy the handle
    checkCudaErrors(cusparseDestroy(handle));

}


int prepare_cusparse_CSR(CSR& cmat, int* csrRowPtr, int *csrColInd, float *csrVal)
{
    if (cmat.fmt != 0)
    {
        std::cout << "ERROR: cusparse_gemm_custom only supports CSR (row-major) " << std::endl;
        return 1;
    }


    //fill csrRowPtr (element i holds number of nonzero entries up to row i)
    csrRowPtr[0] = 0;

    int nnz = 0;
    for (int i = 0; i < cmat.rows; i++)
    {
        nnz += cmat.nzcount[i];
        csrRowPtr[i + 1] = nnz;
    }
    //-------------------------------------------------------------

    //fill csrVal (the nonzero values) and csrColInd (their column indices)
    nnz = 0;
    for (int i = 0; i < cmat.rows; i++)
    {
        std::copy(cmat.ja[i], cmat.ja[i] + cmat.nzcount[i], csrColInd + nnz);
        std::copy(cmat.ma[i], cmat.ma[i] + cmat.nzcount[i], csrVal + nnz);
        nnz += cmat.nzcount[i];

    }
}


/*

//WORK IN PROGRESS
void cublas_blockmat_spsp(const VBS& vbmatA, const VBS& vbmatB, DataT* C, int C_lead_dim, float& dt) {
    //multiplies two VBS matrix (vbmatA and vbmatB); stores into C;
    //vbmatA:       column-major entries storage;
    //              column-major block_storage; 
    //vbmatB:       row-major block_storage; TODO: allow general storage format (implement through cublas transpose)
    //              column-major entries storage;
    //C:            column-major storage; TODO: allow general storage format (implement through cublas transpose)


    cudaDataType_t cuda_type;
    if (typeid(DataT) == typeid(float))    cuda_type = CUDA_R_32F;
    else if (typeid(DataT) == typeid(double))   cuda_type = CUDA_R_64F;
    else if (typeid(DataT) == typeid(int))  cuda_type = CUDA_R_8I;

    cublasGemmAlgo_t cuda_algo = CUBLAS_GEMM_DEFAULT_TENSOR_OP;

    int A_rows = vbmatA.row_part[vbmatA.block_rows];
    int A_cols = vbmatA.col_part[vbmatA.block_cols];

    int B_rows = vbmatB.row_part[vbmatB.block_rows];
    int B_cols = vbmatB.col_part[vbmatB.block_cols];


    int mat_idx = 0; //keeps writing position for mat
    int vbmatA_idx = 0; //keeps reading position for A
    int vbmatB_idx = 0; //keeps reading position for B
    int ja_count_A = 0; //keeps total nonzero blocks count for A;
    int ja_count_B = 0; //keeps total nonzero blocks count for B;

    int C_rows = A_rows;
    int C_cols = B_cols;

    const DataT alpha = 1.0f;
    const DataT beta = 1.0f;

    int tot_nonzero_blocks_A = 0; //index for total nonzero blocks in A 
    int tot_nonzero_blocks_B = 0; //index for total nonzero blocks in B

    int rows_in_block, cols_in_block;
    unsigned int size_block, mem_size_block;

    //TODO: allocate memory on device
    unsigned int size_A = vbmatA.nztot; //total nonzero entries in vbmat
    unsigned int mem_size_A = sizeof(DataT) * size_A;

    unsigned int size_B = vbmatB.nztot;
    unsigned int mem_size_B = sizeof(DataT) * size_B;

    unsigned int size_C = C_rows * C_cols;
    unsigned int mem_size_C = sizeof(DataT) * size_C;

    cublasHandle_t handle;

    checkCudaErrors(cublasCreate(&handle));

    DataT* d_A, * d_B, * d_C;
    checkCudaErrors(cudaMalloc((void**)&d_A, mem_size_A));
    checkCudaErrors(cudaMalloc((void**)&d_B, mem_size_B));
    checkCudaErrors(cudaMalloc((void**)&d_C, mem_size_C));


    //copy to device the vbmat matrix (nonzero blocks are stored consecutively and in column major format)
    checkCudaErrors(cublasSetVector(
        size_A, sizeof(DataT), vbmatA.mab, 1, d_A, 1));

    checkCudaErrors(cublasSetVector(
        size_B, sizeof(DataT), vbmatB.mab, 1, d_B, 1));

    //modify from here

    //initialize cuda events
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);

    //creates streams. Each block row is assigned a different stream.
    cudaStream_t streams[vbmatA.block_rows];
    for (int ib = 0; ib < vbmatA.block_rows; ib++)
    {
        cudaStreamCreate(&(streams[ib]));
    }

    //loop through all blocks
    for (int jb = 0; jb < vbmatA.block_cols; jb++)      //loop horizontally through block columns
    {
        cols_in_block = vbmatA.col_part[jb + 1] - vbmatA.col_part[jb];
        const DataT* d_B_block = d_B + vbmatA.col_part[jb];    //access the block of B that is going to be multiplied with blocks of A in column jb

        for (int nzs = 0; nzs < vbmatA.nzcount[jb]; nzs++)        //loop vertically through nonzero blocks

        {

            int ib = vbmatA.jab[tot_nonzero_blocks];             //the block row position of a nonzero block 
            tot_nonzero_blocks += 1;
            rows_in_block = vbmatA.row_part[ib + 1] - vbmatA.row_part[ib]; //the row height of the block

            cublasSetStream(handle, streams[ib]);               //each block row works on a different stream

            //define the sub-matrices
            const DataT* d_A_block = d_A + vbmat_idx;           //access the block on d_A.
            DataT* d_C_block = d_C + vbmatA.row_part[ib];      //access the block on d_C.




            //multiply the blocks, store result in d_C_block
            checkCudaErrors(
                cublasGemmEx(
                    handle, CUBLAS_OP_N, CUBLAS_OP_N,
                    rows_in_block, B_cols, cols_in_block,           //m, n, k <-- block_A: m*k   block_B: k*n   block_C: m*n
                    &alpha,
                    d_A_block,                                      // blockA device pointer,
                    cuda_type,                                      // blockA datatype
                    rows_in_block,                                  // blockA leading dimension
                    d_B_block,                                      // blockB device pointer
                    cuda_type,                                      // blockB datatype
                    B_rows,                                         // leading dimension
                    &beta,
                    d_C_block, cuda_type,                           // blockC device pointer, blockC type
                    C_rows,
                    cuda_type,                                      // compute_type
                    cuda_algo)
            );
            vbmat_idx += rows_in_block * cols_in_block;
        }

    }

    //record the elapsed time onto dt
//    cudaDeviceSynchronize();
    cudaEventRecord(stop, 0);

    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&dt, start, stop);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);



    //let each stream copy the relevant C block from device
    for (int ib = 0; ib < vbmatA.block_rows; ib++)
    {
        cublasSetStream(handle, streams[ib]);
        rows_in_block = vbmatA.row_part[ib + 1] - vbmatA.row_part[ib];
        checkCudaErrors(cublasGetMatrix(
            rows_in_block, C_cols, sizeof(DataT), d_C + vbmatA.row_part[ib], C_rows, C + vbmatA.row_part[ib], C_lead_dim));

    }

    cudaDeviceSynchronize();

    checkCudaErrors(cudaFree(d_C));
    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_B));

    checkCudaErrors(cublasDestroy(handle));

}


*/