//CUDA Utilities and system includes
#include <assert.h>

// CUDA runtime
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>

#include "driver_types.h"
#include "helper_cuda.h"

//others
#include <stdio.h>
#include <iostream>
#include <typeinfo>

//includes from SPARTA
#include "cuda_utilities.h"
#include "matrices.h"
#include "utilities.h"

/*ONLY WORKS WITH 
    DataT = float, double, int;
    intT = int
*/

void cublas_blockmat_multiply(const VBR& vbmatA, DataT* B, int B_rows, int B_lead_dim, DataT_C* C, int C_lead_dim, float& dt, int n_streams)
{
//multiplies a dense matrix (B) and a VBS matrix (vbmatA); stores B*A into (C)
    //vbmatA:       column-major entries (in-block) storage;
    //              row-major block storage; 
    //B:            column-major storage; TODO: allow general storage format (implement through cublas transpose)
    //C:            column-major storage; TODO: allow general storage format (implement through cublas transpose)


    cudaDataType_t data_type_AB;
    cudaDataType_t data_type_C;
    cublasComputeType_t compute_type;

    if (typeid(DataT) == typeid(int8_t))
    {
        data_type_AB = CUDA_R_8I;
        data_type_C = CUDA_R_32I;
        compute_type = CUBLAS_COMPUTE_32I;
    }
    else if (typeid(DataT) == typeid(float))
    {
        data_type_AB = CUDA_R_32F;
        data_type_C = CUDA_R_32F;
        compute_type = CUBLAS_COMPUTE_32F;
    }
    else
    {
        std::cout << "WARNING! Unsopported multiplication type in cublas_blockmat_multiply(). Check matrices.h" << std::endl;
    }


    cublasGemmAlgo_t cuda_algo = CUBLAS_GEMM_DEFAULT_TENSOR_OP;

    intT A_rows = vbmatA.rows;
    intT A_cols = vbmatA.cols;

    int B_cols = A_rows;
    int C_rows = B_rows;
    int C_cols = A_cols;

    const DataT_C alpha = 1;
    const DataT_C beta = 1;

    //allocate memory on device
    intT size_A = vbmatA.nztot; //total nonzero entries in vbmat
    intT mem_size_A = sizeof(DataT) * size_A;

    intT size_B = B_rows * B_cols;
    intT mem_size_B = sizeof(DataT) * size_B;

    intT size_C = C_rows * C_cols;
    intT mem_size_C = sizeof(DataT_C) * size_C;

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

    
    if (n_streams > vbmatA.block_cols) n_streams = vbmatA.block_cols;
    cudaStream_t streams[n_streams];
    for (intT jb = 0; jb < n_streams; jb++)
    {
        cudaStreamCreate(&(streams[jb]));
    }


    intT mat_idx = 0; //keeps writing position for mat
    intT vbmat_idx = 0; //keeps reading position for vbmat 
    intT ja_count = 0; //keeps total nonzero blocks count;
    intT tot_nonzero_blocks = 0; //index for total nonzero blocks
    intT rows_in_block;
    intT size_block, mem_size_block;
    intT* jab_loc = vbmatA.jab;


    //loop through all blocks
    for(intT ib = 0; ib < vbmatA.block_rows; ib++ )      //loop horizontally through block rows
    {
        rows_in_block = vbmatA.row_part[ib + 1] - vbmatA.row_part[ib]; //the row height of the block
        
        const DataT* d_B_block = d_B + B_rows*vbmatA.row_part[ib];    //access the vertical block of B that is going to be multiplied with blocks of A in block-row ib

        for(intT nzs = 0; nzs < vbmatA.nzcount[ib]; nzs++)        //loop horizontally through nonzero blocks

        {
            intT jb = *jab_loc;             //the block row position of a nonzero block 

            cublasSetStream(handle, streams[jb%n_streams]);               //each stream handles a separate block-column

            //define the sub-matrices
	        const DataT* d_A_block = d_A + vbmat_idx;           //access the block on d_A.
            DataT* d_C_block = d_C + vbmatA.block_col_size*jb*B_rows ;      //access the block on d_C.


            //multiply the blocks, store result in d_C_block
            {
                cublasOperation_t transa = CUBLAS_OP_N, transb = CUBLAS_OP_N;
                int m = B_rows, n = vbmatA.block_col_size, k = rows_in_block;   // TEST k = vbmatA.block_col_size --> B_rows
                int lda = rows_in_block, ldb = B_rows, ldc = C_rows;
                cudaDataType_t Atype = data_type_AB, Btype = data_type_AB, Ctype = data_type_C;

                cublasStatus_t err = cublasGemmEx(
                    handle, CUBLAS_OP_N, CUBLAS_OP_N,
                    B_rows, vbmatA.block_col_size, rows_in_block,           //m, n, k <-- block_B: m*k   block_A: k*n   block_C: m*n
                    &alpha,
                    d_B_block,                                     // blockA device pointer,
                    data_type_AB,                                      // blockA datatype
                    B_rows,                                  // blockA leading dimension
                    d_A_block,                                    // blockB device pointer
                    data_type_AB,                                      // blockB datatype
                    rows_in_block,                                         // leading dimension
                    &beta,
                    d_C_block, Ctype,                           // blockC device pointer, blockC type
                    ldc,
                    compute_type,                                      // compute_type
                    cuda_algo);

                checkCudaErrors( err );
                if (err != CUBLAS_STATUS_SUCCESS) {
                    if (err == CUBLAS_STATUS_INVALID_VALUE) {
                        std::cout << "m = " << m << " n = " << n << " k = " << k << std::endl;
                        std::cout << "lda = " << lda << " >= max(1, " << m << ")? " << std::endl;
                        std::cout << "ldb = " << ldb << " >= max(1, " << k << ")? " << std::endl;
                        std::cout << "ldc = " << ldc << " >= max(1, " << n << ")? " << std::endl;
                    }
                    exit(42);
                }
            }

            
            //move mab and jab pointers forward
            vbmat_idx += rows_in_block*vbmatA.block_col_size;
            jab_loc++;

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
    for (int jb = 0; jb < vbmatA.block_cols; jb++)
    {
        stream_id = jb % n_streams;
        cublasSetStream(handle, streams[stream_id]);
        checkCudaErrors(cublasGetMatrix(
            C_rows, vbmatA.block_col_size, sizeof(DataT_C), d_C + C_rows*jb*vbmatA.block_col_size, C_rows, C + C_rows*jb*vbmatA.block_col_size, C_lead_dim));

    }

    cudaDeviceSynchronize();

    checkCudaErrors(cudaFree(d_C));
    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_B));

    checkCudaErrors(cublasDestroy(handle));
}


void cublas_blockmat_multiplyAB(const VBR& vbmatA, DataT* B, int B_cols, DataT_C* C, float& dt, int n_streams)
{
//multiplies a VBS matrix (vbmatA) and dense matrix (B); stores A*B into (C)
    //vbmatA:       column-major entries (in-block) storage;
    //              row-major block storage; 
    //B:            column-major storage; TODO: allow general storage format (implement through cublas transpose)
    //C:            column-major storage; TODO: allow general storage format (implement through cublas transpose)


    cudaDataType_t data_type_AB;
    cudaDataType_t data_type_C;
    cublasComputeType_t compute_type;

    if (typeid(DataT) == typeid(int8_t))
    {
        data_type_AB = CUDA_R_8I;
        data_type_C = CUDA_R_32I;
        compute_type = CUBLAS_COMPUTE_32I;
    }
    else if (typeid(DataT) == typeid(float))
    {
        data_type_AB = CUDA_R_32F;
        data_type_C = CUDA_R_32F;
        compute_type = CUBLAS_COMPUTE_32F;
    }
    else
    {
        std::cout << "WARNING! Unsopported multiplication type in cublas_blockmat_multiply(). Check matrices.h" << std::endl;
    }


    cublasGemmAlgo_t cuda_algo = CUBLAS_GEMM_DEFAULT_TENSOR_OP;

    intT A_rows = vbmatA.rows;
    intT A_cols = vbmatA.cols;

    int B_rows = A_cols;
    int C_rows = A_rows;
    int C_cols = B_cols;

    const DataT_C alpha = 1;
    const DataT_C beta = 1;

    //allocate memory on device
    intT size_A = vbmatA.nztot; //total nonzero entries in vbmat
    intT mem_size_A = sizeof(DataT) * size_A;

    intT size_B = B_rows * B_cols;
    intT mem_size_B = sizeof(DataT) * size_B;

    intT size_C = C_rows * C_cols;
    intT mem_size_C = sizeof(DataT_C) * size_C;

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
        B_rows, B_cols, sizeof(DataT), B, B_rows, d_B, B_rows));

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


    intT mat_idx = 0; //keeps writing position for mat
    intT vbmat_idx = 0; //keeps reading position for vbmat 
    intT ja_count = 0; //keeps total nonzero blocks count;
    intT tot_nonzero_blocks = 0; //index for total nonzero blocks
    intT rows_in_block;
    intT size_block, mem_size_block;
    intT* jab_loc = vbmatA.jab;

    //loop through all blocks
    for(intT ib = 0; ib < vbmatA.block_rows; ib++ )      //loop horizontally through block rows
    {
        rows_in_block = vbmatA.row_part[ib + 1] - vbmatA.row_part[ib]; //the row height of the block
        

        cublasSetStream(handle, streams[ib%n_streams]);               //each stream handles a separate block-row

        for(intT nzs = 0; nzs < vbmatA.nzcount[ib]; nzs++)        //loop horizontally through nonzero blocks

        {
            intT jb = *jab_loc;             //the block row position of a nonzero block 

            const DataT* d_B_block = d_B + vbmatA.block_col_size*jb;    //access the vertical block of B that is going to be multiplied with blocks of A in block-row ib

            //define the sub-matrices
	        const DataT* d_A_block = d_A + vbmat_idx;           //access the block on d_A.
            DataT* d_C_block = d_C + vbmatA.row_part[ib] ;      //access the block on d_C.


            //multiply the blocks, store result in d_C_block
            checkCudaErrors(
                cublasGemmEx(
                    handle, CUBLAS_OP_N, CUBLAS_OP_N,
                    B_rows, vbmatA.block_col_size, rows_in_block,           //m, n, k <-- block_B: m*k   block_A: k*n   block_C: m*n
                    &alpha,
                    d_A_block,                                     // blockA device pointer,
                    data_type_AB,                                      // blockA datatype
                    rows_in_block,                                  // blockA leading dimension
                    d_B_block,                                    // blockB device pointer
                    data_type_AB,                                      // blockB datatype
                    B_rows,                                       // leading dimension
                    &beta,
                    d_C_block, data_type_C,                           // blockC device pointer, blockC type
                    C_rows,
                    compute_type,                                      // compute_type
                    cuda_algo)
            );                                       
            
            //move mab and jab pointers forward
            vbmat_idx += rows_in_block*vbmatA.block_col_size;
            jab_loc++;

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
            rows_in_block, C_cols, sizeof(DataT_C), d_C + vbmatA.row_part[ib], C_rows, C + vbmatA.row_part[ib], C_rows));

    }

    cudaDeviceSynchronize();

    checkCudaErrors(cudaFree(d_C));
    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_B));

    checkCudaErrors(cublasDestroy(handle));
}




/*
//Matrix-Matrix multiplication with cublas. A,B,C are in column-major order.
//Matrix A and B are in host
//Matrix d_C is in device to allow for accumulation of results
int cublas_gemm_custom(const DataT* A, unsigned int A_rows, unsigned int A_cols, unsigned int lda, const DataT* B, unsigned int B_cols, unsigned int ldb, DataT_C* C, unsigned int ldc, const DataT_C alpha, const DataT_C beta, float& dt)
{
    cudaDataType_t data_type_AB;
    cudaDataType_t data_type_C;
    cublasComputeType_t compute_type;

    if (typeid(DataT) == typeid(int8_t))
    {
        data_type_AB = CUDA_R_8I;
        data_type_C = CUDA_R_32I;
        compute_type = CUBLAS_COMPUTE_32I;
    }
    else if (typeid(DataT) == typeid(float))
    {
        data_type_AB = CUDA_R_32F;
        data_type_C = CUDA_R_32F;
        compute_type = CUBLAS_COMPUTE_32F;
    }
    else
    {
        std::cout << "WARNING! Unsopported multiplication type. Check comp_mats.h" << std::endl;
    }

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
    unsigned int mem_size_C = sizeof(DataT_C) * size_C;




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
            d_A, data_type_AB, A_rows,
            d_B, data_type_AB, B_rows,
            &beta,
            d_C, data_type_C, C_rows,
            compute_type,
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
    checkCudaErrors(cublasGetMatrix(C_rows, C_cols, sizeof(DataT_C), d_C, C_rows, C, C_rows));



    // clean up memory
    checkCudaErrors(cudaFree(d_C));
    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_B));



    // Destroy the handle
    checkCudaErrors(cublasDestroy(handle));
	
    return 0;
}*/


int cusparse_gemm_custom(int rows, int cols, int nnz, int* csrRowPtr, int* csrColInd, DataT* csrVal, DataT* B, int B_cols, int B_lead_dim, DataT_C* C, int C_lead_dim, const DataT_C alpha, const DataT_C beta, float& dt)
{

    cudaDataType_t data_type_AB;
    cudaDataType_t data_type_C;

    if (typeid(DataT) == typeid(int8_t))
    {
        data_type_AB = CUDA_R_8I;
        data_type_C = CUDA_R_32I;
    }
    else if (typeid(DataT) == typeid(float))
    {
        data_type_AB = CUDA_R_32F;
        data_type_C = CUDA_R_32F;
    }
    else
    {
        std::cout << "WARNING! Unsopported multiplication type. Check comp_mats.h" << std::endl;
    }


    //allocate memory on device
    unsigned int mem_size_csrVal = sizeof(DataT) * nnz;
    unsigned int mem_size_csrColInd = sizeof(int) * nnz;
    unsigned int mem_size_csrRowPtr = sizeof(int) * (rows + 1);

    unsigned int B_rows = cols;
    unsigned int size_B = B_rows * B_cols;
    unsigned int mem_size_B = sizeof(DataT) * size_B;
    
    unsigned int C_rows = rows;
    unsigned int C_cols = B_cols;
    unsigned int size_C = C_rows * C_cols;
    unsigned int mem_size_C = sizeof(DataT_C) * size_C;

    // allocate device memory
    int* d_RowPtr, * d_ColInd;
    DataT* d_Val;

    checkCudaErrors(cudaMalloc((void**)&d_RowPtr, mem_size_csrRowPtr));
    checkCudaErrors(cudaMalloc((void**)&d_ColInd, mem_size_csrColInd));
    checkCudaErrors(cudaMalloc((void**)&d_Val, mem_size_csrVal));

    DataT* d_B;
    DataT_C* d_C;
    checkCudaErrors(cudaMalloc((void**)&d_B, mem_size_B));
    checkCudaErrors(cudaMalloc((void**)&d_C, mem_size_C));

    //copy arrays and matrices to device
    // -----------------------
//     checkCudaErrors(cublasSetVector(
//         nnz, sizeof(DataT),
//         csrVal, 1, d_Val, 1));
    // >>>>>>>>>>>>>>>>>>>>>>>
    checkCudaErrors( cudaMemcpy(d_Val, csrVal, mem_size_csrVal, cudaMemcpyHostToDevice) );
    // -----------------------

    checkCudaErrors(cublasSetVector(
        nnz, sizeof(int),
        csrColInd, 1, d_ColInd, 1));

    checkCudaErrors(cublasSetVector(
        (rows + 1), sizeof(int),
        csrRowPtr, 1, d_RowPtr, 1));

    checkCudaErrors(cublasSetMatrix(
        B_rows, B_cols, sizeof(DataT), B, B_lead_dim, d_B, B_rows));


    if (beta != 0)
    {
        checkCudaErrors(cublasSetMatrix(
            C_rows, C_cols, sizeof(DataT_C), C, C_lead_dim, d_C, C_rows));
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
            data_type_AB
            )
        );

    checkCudaErrors(
        cusparseCreateDnMat(&matB, cols, B_cols, B_lead_dim, d_B,
            data_type_AB, CUSPARSE_ORDER_COL));

    checkCudaErrors(
        cusparseCreateDnMat(&matC, rows, B_cols, C_lead_dim, d_C,
            data_type_C, CUSPARSE_ORDER_COL));

    size_t bufferSize;

    checkCudaErrors(cusparseSpMM_bufferSize(
        handle,
        CUSPARSE_OPERATION_NON_TRANSPOSE,
        CUSPARSE_OPERATION_NON_TRANSPOSE,
        (void*)&alpha,
        matA,
        matB,
        (void*)&beta,
        matC,
        data_type_C,
        CUSPARSE_SPMM_ALG_DEFAULT,
        &bufferSize
    ));


    //initialize cuda events
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    
    
    checkCudaErrors(cusparseSpMM(handle,
        CUSPARSE_OPERATION_NON_TRANSPOSE,
        CUSPARSE_OPERATION_NON_TRANSPOSE,
        &alpha, matA, matB, &beta, matC, data_type_C,
        CUSPARSE_SPMM_ALG_DEFAULT, &bufferSize));

    //record the elapsed time onto dt
    cudaDeviceSynchronize();
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&dt, start, stop);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    // destroy matrix/vector descriptors
    checkCudaErrors(cusparseDestroySpMat(matA));
    checkCudaErrors(cusparseDestroyDnMat(matB));
    checkCudaErrors(cusparseDestroyDnMat(matC));

    // copy result from device to host 
    checkCudaErrors(cublasGetMatrix(C_rows, C_cols, sizeof(DataT_C), d_C, C_rows, C, C_rows));

    // clean up memor;y
    checkCudaErrors(cudaFree(d_B));
    checkCudaErrors(cudaFree(d_C));
    checkCudaErrors(cudaFree(d_Val));
    checkCudaErrors(cudaFree(d_RowPtr));
    checkCudaErrors(cudaFree(d_ColInd));

    // Destroy the handle
    checkCudaErrors(cusparseDestroy(handle));

    return 0;
}


int prepare_cusparse_CSR(CSR& cmat, int* csrRowPtr, int* csrColInd, DataT* csrVal)
{
//     if (cmat.fmt != 0)
//     {
//         std::cout << "ERROR: cusparse_gemm_custom only supports CSR (row-major) " << std::endl;
//         return 1;
//     }

    csrRowPtr = (int*) malloc( (cmat.rows+1) * sizeof(int));

    //fill csrRowPtr (element i holds number of nonzero entries up to row i)
    csrRowPtr[0] = 0;

    int nnz = 0;
    for (int i = 0; i < cmat.rows; i++)
    {
        nnz += cmat.nzcount[i];
        csrRowPtr[i + 1] = nnz;
    }
    //-------------------------------------------------------------

    csrColInd = (int*) malloc( nnz * sizeof(int));
    csrVal = (DataT*) malloc( nnz * sizeof(DataT));

    //fill csrVal (the nonzero values) and csrColInd (their column indices)
    nnz = 0;
    for (int i = 0; i < cmat.rows; i++)
    {
        std::copy(cmat.ja[i], cmat.ja[i] + cmat.nzcount[i], csrColInd + nnz);
        std::copy(cmat.ma[i], cmat.ma[i] + cmat.nzcount[i], csrVal + nnz);
        nnz += cmat.nzcount[i];

    }

    return 0;
}
