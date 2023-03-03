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
    // --------------------- TEST Error on cudaMemcpy2D ---------------------
//     checkCudaErrors(cublasSetMatrix(
//         B_rows, B_cols, sizeof(DataT), B, B_lead_dim, d_B, B_rows));
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    checkCudaErrors(cudaMemcpy(d_B, B, B_rows * B_cols * sizeof(DataT), cudaMemcpyHostToDevice));
    // ----------------------------------------------------------------------


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
        // --------------------- TEST Error on cudaMemcpy2D ---------------------
//         checkCudaErrors(cublasGetMatrix(
//             C_rows, vbmatA.block_col_size, sizeof(DataT_C), d_C + C_rows*jb*vbmatA.block_col_size, C_rows, C + C_rows*jb*vbmatA.block_col_size, C_lead_dim));
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        checkCudaErrors(cudaMemcpy(C, d_C, C_rows * C_cols * sizeof(DataT_C), cudaMemcpyDeviceToHost));
        // ----------------------------------------------------------------------
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
    // --------------------- TEST Error on cudaMemcpy2D ---------------------
//     checkCudaErrors(cublasSetMatrix(
//         B_rows, B_cols, sizeof(DataT), B, B_rows, d_B, B_rows));
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    checkCudaErrors(cudaMemcpy(d_B, B, B_rows * B_cols * sizeof(DataT), cudaMemcpyHostToDevice));
    // ----------------------------------------------------------------------

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
        // --------------------- TEST Error on cudaMemcpy2D ---------------------
//         checkCudaErrors(cublasGetMatrix(
//             rows_in_block, C_cols, sizeof(DataT_C), d_C + vbmatA.row_part[ib], C_rows, C + vbmatA.row_part[ib], C_rows));
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        checkCudaErrors(cudaMemcpy(C, d_C, C_rows * C_cols * sizeof(DataT_C), cudaMemcpyDeviceToHost));
        // ----------------------------------------------------------------------

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
    
    long int C_rows = rows;
    long int C_cols = B_cols;
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
    checkCudaErrors(cublasSetVector(
        nnz, sizeof(DataT),
        csrVal, 1, d_Val, 1));
    // >>>>>>>>>>>>>>>>>>>>>>>
//     checkCudaErrors( cudaMemcpy(d_Val, csrVal, mem_size_csrVal, cudaMemcpyHostToDevice) );
    // -----------------------
    void *tmp = malloc(sizeof(DataT) * nnz);
    checkCudaErrors( cudaMemcpy(tmp, d_Val, sizeof(DataT) * nnz, cudaMemcpyDeviceToHost) );
    int out = memcmp(csrVal, tmp, sizeof(DataT) * nnz);
    std::cout << "compare between csrVal and tmp is " << out << "\n";
    free(tmp);

    checkCudaErrors(cublasSetVector(
        nnz, sizeof(int),
        csrColInd, 1, d_ColInd, 1));

    tmp = malloc(sizeof(int) * nnz);
    checkCudaErrors( cudaMemcpy(tmp, d_ColInd, sizeof(int) * nnz, cudaMemcpyDeviceToHost) );
    out = memcmp(csrColInd, tmp, sizeof(int) * nnz);
    std::cout << "compare between csrColInd and tmp is " << out << "\n";
    free(tmp);

    checkCudaErrors(cublasSetVector(
        (rows + 1), sizeof(int),
        csrRowPtr, 1, d_RowPtr, 1));

    tmp = malloc(sizeof(int) * (rows +1) );
    checkCudaErrors( cudaMemcpy(tmp, d_RowPtr, sizeof(int) * (rows +1), cudaMemcpyDeviceToHost) );
    out = memcmp(csrRowPtr, tmp, sizeof(int) * nnz);
    std::cout << "compare between csrRowPtr and tmp is " << out << "\n";
    free(tmp);

    checkCudaErrors(cublasSetMatrix(
        B_rows, B_cols, sizeof(DataT), B, B_lead_dim, d_B, B_rows));

    tmp = malloc(sizeof(DataT) * B_rows * B_cols );
    checkCudaErrors( cudaMemcpy(tmp, d_B, sizeof(DataT) * B_rows * B_cols, cudaMemcpyDeviceToHost) );
    out = memcmp(B, tmp, sizeof(DataT) * B_rows * B_cols);
    std::cout << "compare between B and tmp is " << out << "\n";
    free(tmp);


    if (beta != 0)
    {
        checkCudaErrors(cublasSetMatrix(
            C_rows, C_cols, sizeof(DataT_C), C, C_lead_dim, d_C, C_rows));
    }
    tmp = malloc(sizeof(DataT_C) * C_rows * C_cols );
    checkCudaErrors( cudaMemcpy(tmp, d_C, sizeof(DataT_C) * C_rows * C_cols, cudaMemcpyDeviceToHost) );
    out = memcmp(C, tmp, sizeof(DataT_C) * C_rows * C_cols);
    std::cout << "compare between C and tmp is " << out << "\n";
    free(tmp);

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

    long int chk_rows, chk_cols, chk_ld, chk_nnz;
    void *chk_values, *chk_row, *chk_col;
    cudaDataType chk_datatype;
    cusparseIndexType_t chk_datatypecol, chk_datatyperow;
    cudaDataType_t chk_datatypeval;
    cusparseOrder_t chk_order;
    cusparseIndexBase_t chk_base;
    checkCudaErrors( cusparseCsrGet(matA, &chk_rows, &chk_cols, &chk_nnz, &chk_row, &chk_col, &chk_values, &chk_datatyperow, &chk_datatypecol, &chk_base, &chk_datatypeval) );
    std::cout << "| chk_rows | chk_cols |  chk_nnz  | datatype | chk_order" << std::endl;
    printf("%8ld %8ld %8ld\n", chk_rows, chk_cols, chk_nnz);

    checkCudaErrors(
        cusparseCreateDnMat(&matB, cols, B_cols, B_lead_dim, d_B,
            data_type_AB, CUSPARSE_ORDER_COL));

    checkCudaErrors( cusparseDnMatGet(matB, &chk_rows, &chk_cols, &chk_ld, &chk_values, &chk_datatype, &chk_order) );
    std::cout << "| chk_rows | chk_cols |  chk_ld  | datatype | chk_order" << std::endl;
    printf("%8ld %8ld %8ld\n", chk_rows, chk_cols, chk_ld);

    checkCudaErrors(
        cusparseCreateDnMat(&matC, rows, B_cols, C_lead_dim, d_C,
            data_type_C, CUSPARSE_ORDER_COL));

    checkCudaErrors( cusparseDnMatGet(matC, &chk_rows, &chk_cols, &chk_ld, &chk_values, &chk_datatype, &chk_order) );
    std::cout << "| chk_rows | chk_cols |  chk_ld  | datatype | chk_order" << std::endl;
    printf("%8ld %8ld %8ld\n", chk_rows, chk_cols, chk_ld);

    size_t bufferSize = 0;
    void *dBuffer = NULL;

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
    checkCudaErrors( cudaMalloc(&dBuffer, bufferSize) );


    //initialize cuda events
    cudaEvent_t start, stop;
    checkCudaErrors(cudaEventCreate(&start));
    checkCudaErrors(cudaEventCreate(&stop));
    checkCudaErrors(cudaEventRecord(start, 0));
    

    checkCudaErrors(cusparseSpMM(handle,
        CUSPARSE_OPERATION_NON_TRANSPOSE,
        CUSPARSE_OPERATION_NON_TRANSPOSE,
        &alpha, matA, matB, &beta, matC, data_type_C,
        CUSPARSE_SPMM_ALG_DEFAULT, dBuffer));       // We have a BUG here

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
    // -------------------------------------------------
//     checkCudaErrors(cublasGetMatrix(C_rows, C_cols, sizeof(DataT_C), d_C, C_rows, C, C_rows));
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    checkCudaErrors(cudaMemcpy(C, d_C, C_rows * C_cols * sizeof(DataT_C), cudaMemcpyDeviceToHost));
    // -------------------------------------------------

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

void pico_print_SpMMM(const char* Aname, int An, int Am, int Az, int* Arows, int* Acols, DataT* Avals, const char* Bname, int Bn, int Bm, DataT* B, const char* Cname, long int Cn, long int Cm, DataT_C* C) {

    printf("Dims of %s %s %s:\n", Aname, Bname, Cname);
    printf("       %6s %6s %6s\n", "rows", "cols", "nnz");
    if (Arows != NULL) printf(" %5s %6d %6d %6d\n", Aname, An, Am, Az);
    if (B != NULL)     printf(" %5s %6d %6d %6d\n", Bname, Bn, Bm, 0);
    if (C != NULL)     printf(" %5s %6ld %6ld %6d\n", Cname, Cn, Cm, 0);
    printf("\n");


    if (Arows != NULL) {
        printf("Sparse matrix %s:\n", Aname);
        printf("Rows: ");
        for (int i=0; i<(An+1); i++)
            printf("%3d ", Arows[i]);
        printf("\n");
        printf("Cols: ");
        for (int i=0; i<Az; i++)
            printf("%3d ", Acols[i]);
        printf("\n");
        printf("Vals: ");
        for (int i=0; i<Az; i++)
            printf("%3.2f ", Avals[i]);
        printf("\n\n");
    }

    if (B != NULL) {
        printf("Dense matrix %s:\n", Bname);
        for (int i=0; i<Bn; i++) {
            for (int j=0; j<Bm; j++) {
                printf(" %6.3f ", B[i * Bm + j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    if (C != NULL) {
        printf("Dense matrix %s:\n", Cname);
        for (int i=0; i<Cn; i++) {
            for (int j=0; j<Cm; j++) {
                printf(" %6.3f ", C[i * Cm + j]);
            }
            printf("\n");
        }
        printf("\n");
    }

}

void pico_print_SpMMM(const char* Aname, VBR* A, const char* Bname, int Bn, int Bm, DataT* B, const char* Cname, long int Cn, long int Cm, DataT_C* C) {

    printf("Dims of %s %s %s:\n", Aname, Bname, Cname);
    printf("       %6s %6s %6s\n", "rows", "cols", "block_rows x block_cols");
    if (A != NULL)     printf(" %5s %6ld %6ld      %ldx%ld\n", Aname, A->rows, A->cols, A->block_rows, A->block_cols);
    if (B != NULL)     printf(" %5s %6d %6d      %dx%d\n", Bname, Bn, Bm, 0, 0);
    if (C != NULL)     printf(" %5s %6ld %6ld      %dx%d\n", Cname, Cn, Cm, 0, 0);
    printf("\n");


    if (A != NULL) {
        printf("Sparse matrix %s:\n", Aname);
        printf("nzcount: ");
        for (int i=0; i<((A->rows)/(A->block_rows)); i++)
            printf("%3ld ", A->nzcount[i]);
        printf("\n");
        printf("jab: ");
        int k = 0;
        for (int i=0; i<((A->rows)/(A->block_rows)); i++)
            for (int j=0; j<A->nzcount[i]; j++) {
                printf("%3ld ", A->jab[k]);
                k++;
            }
        printf("\n");
        printf("mab: ");
        for (int i=0; i<A->nztot; i++)
            printf("%3.2f ", A->mab[i]);
        printf("\n");
        printf("row_part: ");
        for (int i=0; i<((A->rows)/(A->block_rows) +1); i++)
            printf("%3.2ld ", A->row_part[i]);
        printf("\n\n");
    }

    if (B != NULL) {
        printf("Dense matrix %s:\n", Bname);
        for (int i=0; i<Bn; i++) {
            for (int j=0; j<Bm; j++) {
                printf(" %6.3f ", B[i * Bm + j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    if (C != NULL) {
        printf("Dense matrix %s:\n", Cname);
        for (int i=0; i<Cn; i++) {
            for (int j=0; j<Cm; j++) {
                printf(" %6.3f ", C[i * Cm + j]);
            }
            printf("\n");
        }
        printf("\n");
    }

}

void pico_print_SpMMM(const char* Aname, int rows, int cols, int ell_blocksize, int ellValue_cols, int ellColumnsInd_rows, int ellColumnsInd_cols, int num_blocks, intT* ellColumnsInd, DataT_C* ellValues, const char* Bname, int Bn, int Bm, DataT* B, const char* Cname, long int Cn, long int Cm, DataT_C* C) {

    printf("Dims of %s %s %s:\n", Aname, Bname, Cname);
    printf("       %6s %6s\n", "rows", "cols");
    if (ellColumnsInd != NULL)     printf(" %5s %6d %6d\n", Aname, rows, cols);
    if (B != NULL)     printf(" %5s %6d %6d\n", Bname, Bn, Bm);
    if (C != NULL)     printf(" %5s %6ld %6ld\n", Cname, Cn, Cm);
    printf("\n");

    printf("rows = %d, cols = %d, ell_blocksize = %d, ellValue_cols = %d, ellColumnsInd_rows = %d, ellColumnsInd_cols = %d, num_blocks = %d\n", rows, cols, ell_blocksize, ellValue_cols, ellColumnsInd_rows, ellColumnsInd_cols, num_blocks);

    if (ellColumnsInd != NULL) {
        printf("Blocked-ellpac matrix %s:\n", Aname);
        printf("\tellColInd matrix:\n");
        for (int i=0; i<ellColumnsInd_rows; i++) {
            printf("\t\t");
            for (int j=0; j< ellColumnsInd_cols; j++)
                printf("%ld ", ellColumnsInd[i*ellColumnsInd_cols + j]);
            printf("\n");
        }

        printf("\tellValue matrix:\n");
        int t = 0;
        for (int i=0; i<ellColumnsInd_rows; i++) {
            for (int j=0; j< ellColumnsInd_cols; j++) {
                if (ellColumnsInd[i*ellColumnsInd_cols + j] != -1) {
                    for(int k=0; k<ell_blocksize; k++) {
                        printf("\t\t");
                        for(int h=0; h<ell_blocksize; h++) {
                            printf("%f ", ellValues[t*ell_blocksize*ell_blocksize +  k*ell_blocksize + h]);
                        }
                        printf("\n");
                    }
                    t++;
                } else {
                    printf("\t\t Padding block\n");
                }
                printf("\t\t--------------------\n");
            }
            printf("\t\t~~~~~~~~~~~~~~~~~~~~\n");
        }

        printf("\trow ellColInd matrix:\n");
        for (int i=0; i<ellColumnsInd_cols*ellColumnsInd_rows; i++)
            printf("%ld ", ellColumnsInd[i]);
        printf("\n");

        printf("\trow ellValue matrix:\n");
        for (int i=0; i<num_blocks*ell_blocksize*ell_blocksize; i++)
            printf("%f ", ellValues[i]);
        printf("\n");


    }

    if (B != NULL) {
        printf("Dense matrix %s:\n", Bname);
        for (int i=0; i<Bn; i++) {
            for (int j=0; j<Bm; j++) {
                printf(" %6.3f ", B[i * Bm + j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    if (C != NULL) {
        printf("Dense matrix %s:\n", Cname);
        for (int i=0; i<Cn; i++) {
            for (int j=0; j<Cm; j++) {
                printf(" %6.3f ", C[i * Cm + j]);
            }
            printf("\n");
        }
        printf("\n");
    }

}

void ck_out (const char* name, bool test) {
    if (test)
        printf("\t%s test: TRUE\n", name);
    else
        printf("\t%s test: FALSE\n", name);
}

void check_csr(cusparseSpMatDescr_t csr, const char* Aname, int An, int Am, int Az, int* Arows, int* Acols, DataT* Avals) {
    printf("Checking csr %s:\n", Aname);

    int64_t chk_rows;
    int64_t chk_cols;
    int64_t chk_nnz;
    void  *chk_csrRowOffsets, *h_chk_csrRowOffsets;
    void  *chk_csrColInd, *h_chk_csrColInd;
    void  *chk_csrValues, *h_chk_csrValues;
    cusparseIndexType_t  chk_csrColIndType;
    cusparseIndexType_t  chk_csrRowOffsetsType;
    cusparseIndexBase_t  chk_idxBase;
    cudaDataType chk_valueType;


    cusparseCsrGet(
        csr,
        &chk_rows,
        &chk_cols,
        &chk_nnz,
        &chk_csrRowOffsets,
        &chk_csrColInd,
        &chk_csrValues,
        &chk_csrRowOffsetsType,
        &chk_csrColIndType,
        &chk_idxBase,
        &chk_valueType
    );

    bool ck_dim = ((An == chk_rows) && (Am == chk_cols) && (Az == chk_nnz));
    ck_out ("dimension", ck_dim);

    h_chk_csrRowOffsets = malloc(sizeof(int) * An);
    cudaMemcpy(h_chk_csrRowOffsets, chk_csrRowOffsets, sizeof(int) * An, cudaMemcpyDeviceToHost);
    bool ck_row = (memcmp(h_chk_csrRowOffsets, Arows, sizeof(int) * An) == 0) ? 1 : 0;
    ck_out ("rows", ck_row);

    h_chk_csrColInd = malloc(sizeof(int) * Az);
    cudaMemcpy(h_chk_csrColInd, chk_csrColInd, sizeof(int) * Az, cudaMemcpyDeviceToHost);
    bool ck_col = (memcmp(h_chk_csrColInd, Acols, sizeof(int) * Az) == 0) ? 1 : 0;
    ck_out ("cols", ck_col);

    h_chk_csrValues = malloc(sizeof(DataT) * Az);
    cudaMemcpy(h_chk_csrValues, chk_csrValues, sizeof(DataT) * Az, cudaMemcpyDeviceToHost);
    bool ck_val = (memcmp(h_chk_csrValues, Avals, sizeof(DataT) * Az) == 0) ? 1 : 0;
    ck_out ("vals", ck_val);

    printf("\n");

// ---------- They are only readed, not copyed ----------
//     cudaFree(chk_csrRowOffsets);
//     cudaFree(chk_csrColInd);
//     cudaFree(chk_csrValues);
    free(h_chk_csrRowOffsets);
    free(h_chk_csrColInd);
    free(h_chk_csrValues);

    return;
}

void check_bell(cusparseSpMatDescr_t bell, const char* Aname, int rows, int cols, int ell_blocksize, int ellValues_cols, int ellColumnsInd_rows, int ellColumnsInd_cols, int num_blocks, intT* ellColumnsInd, DataT_C* ellValues) {
    printf("Checking blocked ellpack %s:\n", Aname);

    int64_t chk_rows;
    int64_t chk_cols;
    int64_t chk_ell_blocksize;
    int64_t chk_ellValues_cols;
    void* chk_ellColumnsInd, *h_chk_ellColumnsInd;
    void* chk_ellValues, *h_chk_ellValues;
    cusparseIndexType_t chk_ellIdxType;
    cusparseIndexBase_t chk_idxBase;
    cudaDataType chk_valueType;

    cusparseBlockedEllGet(bell,
        &chk_rows,
        &chk_cols,
        &chk_ell_blocksize,
        &chk_ellValues_cols,
        &chk_ellColumnsInd,
        &chk_ellValues,
        &chk_ellIdxType,
        &chk_idxBase,
        &chk_valueType);

    bool ck_dim = ((rows == chk_rows) && (cols == chk_cols) && (ell_blocksize == chk_ell_blocksize) && (ellValues_cols == chk_ellValues_cols)), ck_col, ck_val;
    ck_out ("dimension", ck_dim);

    if (ck_dim) {
        h_chk_ellColumnsInd = malloc(sizeof(intT) * (ellColumnsInd_cols*ellColumnsInd_rows));
        cudaMemcpy(h_chk_ellColumnsInd, chk_ellColumnsInd, sizeof(intT) * (ellColumnsInd_cols*ellColumnsInd_rows), cudaMemcpyDeviceToHost);
        ck_col = (memcmp(h_chk_ellColumnsInd, ellColumnsInd, sizeof(intT) * (ellColumnsInd_cols*ellColumnsInd_rows)) == 0) ? 1 : 0;
        ck_out ("ellColumnsInd_cols", ck_col);
    } else {
        printf("(rows == chk_rows): (%d, %ld) && (cols == chk_cols): (%d, %ld) && (ell_blocksize == chk_ell_blocksize): (%d, %ld) && (ellValues_cols == chk_ellValues_cols): (%d, %ld)\n", rows, chk_rows, cols, chk_cols, ell_blocksize, chk_ell_blocksize, ellValues_cols, chk_ellValues_cols);
    }

    if (ck_dim) {
        h_chk_ellValues = malloc(sizeof(DataT_C) * (num_blocks*ell_blocksize*ell_blocksize));
        cudaMemcpy(h_chk_ellValues, chk_ellValues, sizeof(DataT_C) * (num_blocks*ell_blocksize*ell_blocksize), cudaMemcpyDeviceToHost);
        ck_val = (memcmp(h_chk_ellValues, ellValues, sizeof(DataT_C) * (num_blocks*ell_blocksize*ell_blocksize)) == 0) ? 1 : 0;
        ck_out ("ell_val", ck_val);
    }

    if (ck_col == false || ck_val == false) {
        printf("h_chk_matrix:\n");
        pico_print_SpMMM("BEL_A", chk_rows, chk_cols, chk_ell_blocksize, chk_ellValues_cols, ellColumnsInd_rows, ellColumnsInd_cols, num_blocks, (intT*)h_chk_ellColumnsInd, (DataT_C*)h_chk_ellValues, "NULL", 0, 0, NULL, "NULL", 0, 0, NULL);
    }

    printf("\n");

// ---------- They are only readed, not copyed ----------
//     cudaFree(chk_ellColumnsInd);
//     cudaFree(chk_ellValues);
    free(h_chk_ellColumnsInd);
    free(h_chk_ellValues);

    return;
}

void check_dnmat(cusparseDnMatDescr_t dnmat, const char* Bname, int Bn, int Bm, DataT* B) {
    printf("Checking dense matrix %s:\n", Bname);


    int64_t chk_rows;
    int64_t chk_cols;
    int64_t chk_ld;
    void  *chk_vals, *h_chk_vals;
    cudaDataType chk_valueType;
    cusparseOrder_t chk_order;

    cusparseDnMatGet(dnmat,
        &chk_rows,
        &chk_cols,
        &chk_ld,
        &chk_vals,
        &chk_valueType,
        &chk_order
    );



    bool ck_dim = ((Bn == chk_rows) && (Bm == chk_cols) && (Bm == chk_ld));
    ck_out ("dimension", ck_dim);

    h_chk_vals = malloc(sizeof(DataT) * Bn * Bm);
    cudaMemcpy(h_chk_vals, chk_vals, sizeof(DataT) * Bn * Bm, cudaMemcpyDeviceToHost);
    bool ck_val = (memcmp(h_chk_vals, B, sizeof(DataT) * Bn * Bm) == 0) ? 1 : 0;
    ck_out ("vals", ck_val);

    printf("\n");


// ---------- They are only readed, not copyed ----------
//     cudaFree(chk_vals);
    free(h_chk_vals);

    return;
}

// #define BYPASS
#ifndef BYPASS
int cusparse_gemm_custom2(int rows, int cols, int nnz, int* csrRowPtr, int* csrColInd, DataT* csrVal, DataT* B, int B_cols, int B_lead_dim, DataT_C* C, int C_lead_dim, const DataT_C alpha, const DataT_C beta, float& dt)
#else
int cusparse_gemm_custom2(int rows, int cols, int nnz, int* csrRowPtr0, int* csrColInd0, DataT* csrVal0, DataT* B0, int B_cols, int B_lead_dim, DataT_C* C0, int C_lead_dim, const DataT_C alpha, const DataT_C beta, float& dt)
#endif
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

    long int C_rows = rows;
    long int C_cols = B_cols;
    unsigned int size_C = C_rows * C_cols;
    unsigned int mem_size_C = sizeof(DataT_C) * size_C;

    // allocate device memory
    int *d_RowPtr, *d_ColInd;
    DataT *d_Val;

#ifdef BYPASS

    pico_print_SpMMM("A0", rows, cols, nnz, csrRowPtr0, csrColInd0, csrVal0, "B0", B_rows, B_cols, B0, "C0", C_rows, C_cols, C0);

    /*int*/   rows      = 4;
    /*int*/   cols      = 4;
    /*int*/   nnz       = 9;
    /*unsigned int*/   B_rows      = cols;
    /*int*/   B_cols      = 3;
    /*int*/   B_lead_dim          = B_cols;
    /*int*/   C_lead_dim          = B_cols;
    /*int*/   mem_size_B          = B_rows * B_cols * sizeof(DataT);
    /*int*/   mem_size_C          = C_rows * B_cols * sizeof(DataT_C);

// ------------------------ Test ----------------------------
    int   csrRowPtr[] = { 0, 3, 4, 7, 9 };
    int   csrColInd[]    = { 0, 2, 3, 1, 0, 2, 3, 1, 3 };
    DataT csrVal[]     = { 1.0f, 2.0f, 3.0f, 4.0f, 5.0f,
                              6.0f, 7.0f, 8.0f, 9.0f };
    DataT B[]            = { 1.0f,  5.0f,  9.0f,
                             2.0f,  6.0f, 10.0f,
                             3.0f,  7.0f, 11.0f,
                             4.0f,  8.0f, 12.0f };
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//     nnz = 12;
//     /*int*/   mem_size_B          = B_rows * B_cols * sizeof(DataT);
//     /*int*/   mem_size_C          = C_rows * C_cols * sizeof(DataT_C);
//
//     int   csrRowPtr[] = { 0, 0, 4, 8, 12 };
//     int   csrColInd[]    = { 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3 };
//     DataT csrVal[]     = { 1.0f, 1.0f, 1.0f, 1.0f,
//                            2.0f, 2.0f, 2.0f, 2.0f,
//                            3.0f, 3.0f, 3.0f, 3.0f };
//
//     DataT B[]            = { 1.0f,  0.1f,  0.01f,
//                              1.0f,  0.1f,  0.01f,
//                              1.0f,  0.1f,  0.01f,
//                              1.0f,  0.1f,  0.01f };
// ----------------------------------------------------------

    DataT_C C[]            = { 0.0f, 0.0f, 0.0f,
                               0.0f, 0.0f, 0.0f,
                               0.0f, 0.0f, 0.0f,
                               0.0f, 0.0f, 0.0f };
    DataT_C hC_result[]     = { 19.0f,  8.0f,  51.0f,  52.0f,
                              43.0f, 24.0f, 123.0f, 120.0f,
                              67.0f, 40.0f, 195.0f, 188.0f };
//     float alpha           = 1.0f;
//     float beta            = 0.0f;

    // -----------------
    C_rows = rows;
    C_cols = B_cols;

    mem_size_csrVal = sizeof(DataT) * nnz;
    mem_size_csrColInd = sizeof(int) * nnz;
    mem_size_csrRowPtr = sizeof(int) * (rows + 1);

#endif

#ifdef PICO_DEBUG
    printf("alpha = %f, beta = %f\n", alpha, beta);
    pico_print_SpMMM("A", rows, cols, nnz, csrRowPtr, csrColInd, csrVal, "B", B_rows, B_cols, B, "C", C_rows, C_cols, C);
#endif

    checkCudaErrors(cudaMalloc((void**)&d_RowPtr, mem_size_csrRowPtr));
    checkCudaErrors(cudaMalloc((void**)&d_ColInd, mem_size_csrColInd));
    checkCudaErrors(cudaMalloc((void**)&d_Val, mem_size_csrVal));

    DataT* d_B;
    DataT_C* d_C;
    checkCudaErrors(cudaMalloc((void**)&d_B, mem_size_B));
    checkCudaErrors(cudaMalloc((void**)&d_C, mem_size_C));

    //copy arrays and matrices to device
    checkCudaErrors( cudaMemcpy(d_Val, csrVal, mem_size_csrVal, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(d_ColInd, csrColInd, mem_size_csrColInd, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(d_RowPtr, csrRowPtr, mem_size_csrRowPtr, cudaMemcpyHostToDevice) );

    checkCudaErrors( cudaMemcpy(d_B, B, mem_size_B, cudaMemcpyHostToDevice) );

    if (beta != 0)
    {
        checkCudaErrors( cudaMemcpy(d_C, C, mem_size_C, cudaMemcpyHostToDevice) );
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

#ifdef PICO_DEBUG
    check_csr(matA, "A", rows, cols, nnz, csrRowPtr, csrColInd, csrVal);
#endif

    checkCudaErrors(
        cusparseCreateDnMat(&matB, cols, B_cols, B_lead_dim, d_B,
            data_type_AB, CUSPARSE_ORDER_ROW));

#ifdef PICO_DEBUG
    check_dnmat(matB, "matB", B_rows, B_cols, B);
#endif

    checkCudaErrors(
        cusparseCreateDnMat(&matC, rows, B_cols, C_lead_dim, d_C,
            data_type_C, CUSPARSE_ORDER_ROW));

#ifdef PICO_DEBUG
    check_dnmat(matC, "matC", C_rows, C_cols, C);
#endif



    size_t bufferSize = 0;
    void *dBuffer = NULL;

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
    checkCudaErrors( cudaMalloc(&dBuffer, bufferSize) );


    //initialize cuda events
    cudaEvent_t start, stop;
    checkCudaErrors(cudaEventCreate(&start));
    checkCudaErrors(cudaEventCreate(&stop));
    checkCudaErrors(cudaEventRecord(start, 0));


    checkCudaErrors(cusparseSpMM(handle,
        CUSPARSE_OPERATION_NON_TRANSPOSE,
        CUSPARSE_OPERATION_NON_TRANSPOSE,
        &alpha, matA, matB, &beta, matC, data_type_C,
        CUSPARSE_SPMM_ALG_DEFAULT, dBuffer));       // We have a BUG here

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
    // -------------------------------------------------
//     checkCudaErrors(cublasGetMatrix(C_rows, C_cols, sizeof(DataT_C), d_C, C_rows, C, C_rows));
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    checkCudaErrors(cudaMemcpy(C, d_C, mem_size_C, cudaMemcpyDeviceToHost));
    // -------------------------------------------------

#ifdef PICO_DEBUG
    printf("Printing results of %s:\n", __func__);
    pico_print_SpMMM("A", rows, cols, nnz, csrRowPtr, csrColInd, csrVal, "B", B_rows, B_cols, B, "C", C_rows, C_cols, C);
#endif

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

int prepare_cusparse_CSR(CSR& cmat, int** csrRowPtr, int** csrColInd, DataT** csrVal)
{
//     if (cmat.fmt != 0)
//     {
//         std::cout << "ERROR: cusparse_gemm_custom only supports CSR (row-major) " << std::endl;
//         return 1;
//     }

#ifdef PICO_DEBUG
    printf("File %s, line %d: print cmat\n", __FILE__, __LINE__);
    cmat.print(2);
#endif

    *csrRowPtr = (int*) malloc( (cmat.rows+1) * sizeof(int));

    //fill csrRowPtr (element i holds number of nonzero entries up to row i)
    (*csrRowPtr)[0] = 0;

    int nnz = 0;
    for (int i = 0; i < cmat.rows; i++)
    {
        nnz += cmat.nzcount[i];
        (*csrRowPtr)[i + 1] = nnz;
    }
    //-------------------------------------------------------------

    *csrColInd = (int*) malloc( nnz * sizeof(int));
    *csrVal = (DataT*) malloc( nnz * sizeof(DataT));

    //fill csrVal (the nonzero values) and csrColInd (their column indices)
    nnz = 0;
    for (int i = 0; i < cmat.rows; i++)
    {
        std::copy(cmat.ja[i], cmat.ja[i] + cmat.nzcount[i], (*csrColInd) + nnz);
        std::copy(cmat.ma[i], cmat.ma[i] + cmat.nzcount[i], (*csrVal) + nnz);
        nnz += cmat.nzcount[i];

    }

    return 0;
}



int cusparse_gemm_custom_ellpack(int rows, int cols, int A_ell_blocksize, int A_ellValues_cols, int A_ellColInd_cols, int A_ellColInd_rows, int A_num_blocks, intT* A_ellColInd, DataT_C* A_ellValues, DataT* B, int B_cols, int B_lead_dim, DataT_C* C, int C_lead_dim, const DataT_C alpha, const DataT_C beta, float& dt)
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


    unsigned int B_rows = cols;
    unsigned int size_B = B_rows * B_cols;
    unsigned int mem_size_B = sizeof(DataT) * size_B;

// TO ADAPT

    long int C_rows = rows;
    long int C_cols = B_cols;
    unsigned int size_C = C_rows * C_cols;
    unsigned int mem_size_C = sizeof(DataT_C) * size_C;

    // allocate device memory
    intT    *dA_ellColInd;
    DataT_C *dA_ellValues;


#ifdef PICO_DEBUG
    printf("alpha = %f, beta = %f\n", alpha, beta);
//     pico_print_SpMMM("A", rows, cols, nnz, csrRowPtr, csrColInd, csrVal, "B", B_rows, B_cols, B, "C", C_rows, C_cols, C);
#endif

    checkCudaErrors( cudaMalloc((void**) &dA_ellColInd, A_ellColInd_cols * A_ellColInd_rows * sizeof(intT)) );
    checkCudaErrors( cudaMalloc((void**) &dA_ellValues, A_ell_blocksize * A_ell_blocksize * A_num_blocks * sizeof(DataT_C)) );

    DataT* d_B;
    DataT_C* d_C;
    checkCudaErrors(cudaMalloc((void**)&d_B, mem_size_B));
    checkCudaErrors(cudaMalloc((void**)&d_C, mem_size_C));

    //copy arrays and matrices to device
    checkCudaErrors( cudaMemcpy(dA_ellColInd, A_ellColInd, A_ellColInd_rows * A_ellColInd_cols * sizeof(intT), cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(dA_ellValues, A_ellValues, A_num_blocks * A_ell_blocksize * A_ell_blocksize * sizeof(DataT_C), cudaMemcpyHostToDevice) );

    checkCudaErrors( cudaMemcpy(d_B, B, mem_size_B, cudaMemcpyHostToDevice) );

    if (beta != 0)
    {
        checkCudaErrors( cudaMemcpy(d_C, C, mem_size_C, cudaMemcpyHostToDevice) );
    }

    cusparseHandle_t handle;
    cusparseSpMatDescr_t matA;
    cusparseDnMatDescr_t matB, matC;

    checkCudaErrors(cusparseCreate(&handle));

    checkCudaErrors( cusparseCreateBlockedEll(
                                      &matA,
                                      rows, cols, A_ell_blocksize,
                                      A_ellValues_cols, dA_ellColInd, dA_ellValues,
                                      CUSPARSE_INDEX_32I,
                                      CUSPARSE_INDEX_BASE_ZERO, data_type_AB) );

#ifdef PICO_DEBUG
    check_bell(matA, "A", rows, cols, A_ell_blocksize, A_ellValues_cols, A_ellColInd_rows, A_ellColInd_cols, A_num_blocks, A_ellColInd, A_ellValues);
#endif

    checkCudaErrors(
        cusparseCreateDnMat(&matB, cols, B_cols, B_lead_dim, d_B,
            data_type_AB, CUSPARSE_ORDER_ROW));

#ifdef PICO_DEBUG
    check_dnmat(matB, "matB", B_rows, B_cols, B);
#endif

    checkCudaErrors(
        cusparseCreateDnMat(&matC, rows, B_cols, C_lead_dim, d_C,
            data_type_C, CUSPARSE_ORDER_ROW));

#ifdef PICO_DEBUG
    check_dnmat(matC, "matC", C_rows, C_cols, C);
#endif



    size_t bufferSize = 0;
    void *dBuffer = NULL;

    checkCudaErrors( cusparseSpMM_bufferSize(
                                 handle,
                                 CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &alpha, matA, matB, &beta, matC, data_type_C,
                                 CUSPARSE_SPMM_ALG_DEFAULT, &bufferSize) );
    checkCudaErrors( cudaMalloc(&dBuffer, bufferSize) );


    //initialize cuda events
    cudaEvent_t start, stop;
    checkCudaErrors(cudaEventCreate(&start));
    checkCudaErrors(cudaEventCreate(&stop));
    checkCudaErrors(cudaEventRecord(start, 0));

    checkCudaErrors( cusparseSpMM(handle,
                                CUSPARSE_OPERATION_NON_TRANSPOSE,
                                CUSPARSE_OPERATION_NON_TRANSPOSE,
                                &alpha, matA, matB, &beta, matC, data_type_C,
                                CUSPARSE_SPMM_ALG_DEFAULT, dBuffer) );

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
    // -------------------------------------------------
//     checkCudaErrors(cublasGetMatrix(C_rows, C_cols, sizeof(DataT_C), d_C, C_rows, C, C_rows));
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    checkCudaErrors(cudaMemcpy(C, d_C, mem_size_C, cudaMemcpyDeviceToHost));
    // -------------------------------------------------

#ifdef PICO_DEBUG
    printf("Printing results of %s:\n", __func__);
    pico_print_SpMMM("NULL", 0, 0, 0, NULL, NULL, NULL, "B", B_rows, B_cols, B, "C", C_rows, C_cols, C);
#endif

    // clean up memor;y
    checkCudaErrors(cudaFree(d_B));
    checkCudaErrors(cudaFree(d_C));
    checkCudaErrors(cudaFree(dA_ellValues));
    checkCudaErrors(cudaFree(dA_ellColInd));

    // Destroy the handle
    checkCudaErrors(cusparseDestroy(handle));

    return 0;
}

int prepare_cusparse_BLOCKEDELLPACK(VBR *A, int *ell_blocksize, int *ellValue_cols, int *ellColInd_rows, int *ellColInd_cols, int *num_blocks, intT** ellColInd, DataT_C** ellValues)
{
//     if (cmat.fmt != 0)
//     {
//         std::cout << "ERROR: cusparse_gemm_custom only supports CSR (row-major) " << std::endl;
//         return 1;
//     }

#ifdef PICO_DEBUG
    printf("File %s, line %d: print cmat\n", __FILE__, __LINE__);
    A->print(2);
#endif

    *ell_blocksize = A->block_col_size;
    *ellColInd_rows = ((A->rows)/(*ell_blocksize));
    *ellColInd_cols = 0;
    *num_blocks = 0;
    for (int i=0; i<(*ellColInd_rows); i++) {
        *num_blocks += A->nzcount[i];
        if (A->nzcount[i] > (*ellColInd_cols))
            (*ellColInd_cols) = A->nzcount[i];
    }
    *ellValue_cols = (*ellColInd_cols) * (*ell_blocksize);

    *ellColInd = (intT*)malloc(sizeof(intT)*(*ellColInd_cols)*(*ellColInd_rows));
    *ellValues = (DataT_C*)malloc(sizeof(DataT)*((*num_blocks))*(*ell_blocksize)*(*ell_blocksize));

    int i, j, k_col = 0, k_val_i, k_val_j;
    for (i=0; i<(*ellColInd_rows); i++)
        for (j=0; j<(*ellColInd_cols); j++)
            if (j < A->nzcount[i]) {
                (*ellColInd)[i*(*ellColInd_cols)+j] = A->jab[k_col];

                // blocks form colum-major (A->mab) to row-major (*ellValues)
                for (k_val_i=0; k_val_i<(*ell_blocksize); k_val_i++)
                    for (k_val_j=0; k_val_j<(*ell_blocksize); k_val_j++) {
                        (*ellValues)[k_col * (*ell_blocksize)*(*ell_blocksize) + k_val_i * (*ell_blocksize) + k_val_j] =
                            A->mab[k_col * (*ell_blocksize)*(*ell_blocksize) + k_val_j * (*ell_blocksize) + k_val_i];
                    }

                k_col ++;
            } else {
                (*ellColInd)[i*(*ellColInd_cols)+j] = -1; // '-1' means thet the algorithm autmaticly pad it with a compleate zero block
            }

    return 0;
}
