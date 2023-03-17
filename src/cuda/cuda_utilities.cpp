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

/*
void cublas_gemm(DataT *A, intT A_rows, intT A_cols, DataT *B, intT B_cols, float &dt)
{

    //define data types
    cudaDataType_t data_type_AB;
    cudaDataType_t data_type_C;
    cublasComputeType_t compute_type;
    cublasGemmAlgo_t cuda_algo = CUBLAS_GEMM_DEFAULT_TENSOR_OP;

    int B_rows = A_cols;
    int C_rows = A_rows;
    int C_cols = B_cols;

    const DataT_C alpha = 1;
    const DataT_C beta = 1;

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

    //m, n, k <-- block_A: m*k   block_B: k*n   block_C: m*n
    cublasOperation_t transa = CUBLAS_OP_N, transb = CUBLAS_OP_N;
    int m = A_rows, n = B_cols, k = A_cols;
    int lda = A_rows, ldb = B_rows, ldc = C_rows;
    cudaDataType_t Atype = data_type_AB, Btype = data_type_AB, Ctype = data_type_C;

    //allocate memory on device
    intT size_A = A_rows * A_cols; //total nonzero entries in A
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
        size_A, sizeof(DataT), A, 1, d_A, 1));

    checkCudaErrors(cublasSetVector(
        size_B, sizeof(DataT), B, 1, d_B, 1));

    //initialize cuda events
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);

    cublasStatus_t err = cublasGemmEx(
        handle, transa, transb,
        A_rows, B_cols, A_cols,           //m, n, k <-- block_B: m*k   block_A: k*n   block_C: m*n
        &alpha,
        d_A,                                     // blockA device pointer,
        data_type_AB,                                      // blockA datatype
        lda,                                  // blockA leading dimension
        d_B,                                    // blockB device pointer
        data_type_AB,                                      // blockB datatype
        ldb,                                         // leading dimension
        &beta,
        d_C, Ctype,                           // blockC device pointer, blockC type
        ldc,
        compute_type,                                      // compute_type
        cuda_algo);

    checkCudaErrors( err );
    if (err != CUBLAS_STATUS_SUCCESS) 
    {
        if (err == CUBLAS_STATUS_INVALID_VALUE) 
        {
            std::cout << "m = " << m << " n = " << n << " k = " << k << std::endl;
            std::cout << "lda = " << lda << " >= max(1, " << m << ")? " << std::endl;
            std::cout << "ldb = " << ldb << " >= max(1, " << k << ")? " << std::endl;
            std::cout << "ldc = " << ldc << " >= max(1, " << n << ")? " << std::endl;
        }
        exit(42);
    }

    
void cublas_gemm(DataT *A, intT A_rows, intT A_cols, DataT *B, intT B_cols, float &dt)
{

    //define data types
    cudaDataType_t data_type_AB;
    cudaDataType_t data_type_C;
    cublasComputeType_t compute_type;
    cublasGemmAlgo_t cuda_algo = CUBLAS_GEMM_DEFAULT_TENSOR_OP;

    int B_rows = A_cols;
    int C_rows = A_rows;
    int C_cols = B_cols;

    const DataT_C alpha = 1;
    const DataT_C beta = 1;

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

    //m, n, k <-- block_A: m*k   block_B: k*n   block_C: m*n
    cublasOperation_t transa = CUBLAS_OP_N, transb = CUBLAS_OP_N;
    int m = A_rows, n = B_cols, k = A_cols;
    int lda = A_rows, ldb = B_rows, ldc = C_rows;
    cudaDataType_t Atype = data_type_AB, Btype = data_type_AB, Ctype = data_type_C;

    //allocate memory on device
    intT size_A = A_rows * A_cols; //total nonzero entries in A
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
        size_A, sizeof(DataT), A, 1, d_A, 1));

    checkCudaErrors(cublasSetVector(
        size_B, sizeof(DataT), B, 1, d_B, 1));

    //initialize cuda events
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);

    cublasStatus_t err = cublasGemmEx(
        handle, transa, transb,
        m, n, k,           //m, n, k <-- block_B: m*k   block_A: k*n   block_C: m*n
        &alpha,
        d_A,                                     // blockA device pointer,
        data_type_AB,                                      // blockA datatype
        lda,                                  // blockA leading dimension
        d_B,                                    // blockB device pointer
        data_type_AB,                                      // blockB datatype
        ldb,                                         // leading dimension
        &beta,
        d_C, Ctype,                           // blockC device pointer, blockC type
        ldc,
        compute_type,                                      // compute_type
        cuda_algo);

    checkCudaErrors( err );
    if (err != CUBLAS_STATUS_SUCCESS) 
    {
        if (err == CUBLAS_STATUS_INVALID_VALUE) 
        {
            std::cout << "m = " << m << " n = " << n << " k = " << k << std::endl;
            std::cout << "lda = " << lda << " >= max(1, " << m << ")? " << std::endl;
            std::cout << "ldb = " << ldb << " >= max(1, " << k << ")? " << std::endl;
            std::cout << "ldc = " << ldc << " >= max(1, " << n << ")? " << std::endl;
        }
        exit(42);
    }

    //record the elapsed time onto dt
    cudaDeviceSynchronize();
    cudaEventRecord(stop, 0);

    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&dt, start, stop);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    checkCudaErrors(cudaMemcpy(C, d_C, C_rows * C_cols * sizeof(DataT_C), cudaMemcpyDeviceToHost));

    cudaDeviceSynchronize();

    checkCudaErrors(cudaFree(d_C));
    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_B));

    checkCudaErrors(cublasDestroy(handle));


}

}
*/


void cublas_fixed_blocks_multiply(const VBR& vbmatA, DataT* B, int B_cols, DataT_C* C, float& dt, int n_streams)
{
//multiplies a fixed size VBR matrix (vbmatA) and dense matrix (B); stores A*B into (C)
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

    intT block_size = vbmatA.block_col_size;

    intT A_rows = vbmatA.rows;
    intT A_cols = vbmatA.cols;

    intT B_rows = A_cols;
    intT C_rows = A_rows;
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

    intT tot_nz_blocks = std::accumulate(vbmatA.nzcount, vbmatA.nzcount + vbmatA.block_rows);
    intT max_blocks_per_stream = tot_nz_blocks/n_streams + 1;


    intT mat_idx = 0; //keeps writing position for mat
    intT vbmat_idx = 0; //keeps reading position for vbmat 
    intT ja_count = 0; //keeps total nonzero blocks count;
    intT tot_nonzero_blocks = 0; //index for total nonzero blocks
    intT rows_in_block;
    intT nz_blocks = 0;
    intT size_block = block_size*block_size;
    intT* jab_loc = vbmatA.jab;


    intT current_stream = 0;
    intT blocks_in_stream = 0;
    cublasSetStream(handle, streams[0]);               //each stream handles at most max_blocks_per_stream of block_rows

    //loop through all blocks
    for(intT ib = 0; ib < vbmatA.block_rows; ib++ )      //loop horizontally through block rows
    {
        rows_in_block = block_size; //the row height of the block
        nz_blocks = vbmatA.nzcount[ib];
        
        if (blocks_in_stream + nz_blocks > max_blocks_per_stream)
        {
            current_stream++;
            if(current_stream < n_streams) blocks_in_stream = nz_blocks;              //only reset the counter if we did not circle back to the first stream
            cublasSetStream(handle, streams[current_stream%n_streams]);               //each stream handles a separate block-row
        }
        else
        {
            blocks_in_stream += nz_blocks;
        }


        for(intT nzs = 0; nzs < vbmatA.nzcount[ib]; nzs++)        //loop horizontally through nonzero blocks
        {
            intT jb = *jab_loc;             //the block row position of a nonzero block 

            const DataT* d_B_block = d_B + vbmatA.block_col_size*jb;    //access the vertical block of B that is going to be multiplied with blocks of A in block-row ib

            //define the sub-matrices
	        const DataT* d_A_block = d_A + vbmat_idx;           //access the block on d_A.
            DataT* d_C_block = d_C + vbmatA.row_part[ib] ;      //access the block on d_C.

            int m = rows_in_block, n = B_cols, k = vbmatA.block_col_size;
            int lda = rows_in_block, ldb = B_rows, ldc = C_rows;

            //multiply the blocks, store result in d_C_block
            checkCudaErrors(
                cublasGemmEx(
                    handle, CUBLAS_OP_N, CUBLAS_OP_N,
                    m,n,k,                                               //m, n, k <-- block_A: m*k   block_B: k*n   block_C: m*n
                    &alpha,
                    d_A_block,                                          // blockA device pointer,
                    data_type_AB,                                       // blockA datatype
                    lda,                                      // blockA leading dimension
                    d_B_block,                                          // blockB device pointer
                    data_type_AB,                                       // blockB datatype
                    ldb,                                             // B leading dimension
                    &beta,
                    d_C_block, data_type_C,                             // blockC device pointer, blockC type
                    ldc,                                             // C leading dimension
                    compute_type,                                       // compute_type
                    cuda_algo)
            );                                       
            
            //move mab and jab pointers forward
            vbmat_idx += size_block;
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
        checkCudaErrors(cublasGetMatrixAsync(
            rows_in_block, C_cols, sizeof(DataT_C), d_C + vbmatA.row_part[ib], C_rows, C + vbmatA.row_part[ib], C_rows, streams[stream_id]));
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

    intT B_rows = A_cols;
    intT C_rows = A_rows;
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
        DataT* d_C_block = d_C + vbmatA.row_part[ib] ;      //access the block on d_C.
       

        cublasSetStream(handle, streams[ib%n_streams]);               //each stream handles a separate block-row

        for(intT nzs = 0; nzs < vbmatA.nzcount[ib]; nzs++)        //loop horizontally through nonzero blocks

        {
            intT jb = *jab_loc;             //the block row position of a nonzero block 

            const DataT* d_B_block = d_B + vbmatA.block_col_size*jb;    //access the vertical block of B that is going to be multiplied with blocks of A in block-row ib

            //define the sub-matrices
	        const DataT* d_A_block = d_A + vbmat_idx;           //access the block on d_A.

            int m = rows_in_block, n = B_cols, k = vbmatA.block_col_size;
            int lda = rows_in_block, ldb = B_rows, ldc = C_rows;

            //multiply the blocks, store result in d_C_block
            checkCudaErrors(
                cublasGemmEx(
                    handle, CUBLAS_OP_N, CUBLAS_OP_N,
                    m,n,k,                                               //m, n, k <-- block_A: m*k   block_B: k*n   block_C: m*n
                    &alpha,
                    d_A_block,                                          // blockA device pointer,
                    data_type_AB,                                       // blockA datatype
                    lda,                                      // blockA leading dimension
                    d_B_block,                                          // blockB device pointer
                    data_type_AB,                                       // blockB datatype
                    ldb,                                             // B leading dimension
                    &beta,
                    d_C_block, data_type_C,                             // blockC device pointer, blockC type
                    ldc,                                             // C leading dimension
                    compute_type,                                       // compute_type
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
        checkCudaErrors(cublasGetMatrixAsync(
            rows_in_block, C_cols, sizeof(DataT_C), d_C + vbmatA.row_part[ib], C_rows, C + vbmatA.row_part[ib], C_rows, streams[stream_id]));
    }

    cudaDeviceSynchronize();

    checkCudaErrors(cudaFree(d_C));
    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_B));

    checkCudaErrors(cublasDestroy(handle));
}



void cublas_blockmat_multiplyBA(const VBR& vbmatA, DataT* B, int B_rows, DataT_C* C, float& dt, int n_streams)
{
//multiplies a VBS matrix (vbmatA) and dense matrix (B); stores B*A into (C)
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
    intT B_cols = A_rows;
    intT C_rows = B_rows;
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
        
        const DataT* d_B_block = d_B + vbmatA.block_col_size*ib;    //access the vertical block of B that is going to be multiplied with blocks of A in block-row ib

        for(intT nzs = 0; nzs < vbmatA.nzcount[ib]; nzs++)        //loop horizontally through nonzero blocks
        {
            DataT* d_C_block = d_C + vbmatA.block_col_size*C_rows;      //access the block on d_C.
            
            intT jb = *jab_loc;             //the block row position of a nonzero block 

            cublasSetStream(handle, streams[jb%n_streams]);               //each stream handles a separate block-row
            
            //define the sub-matrices
	        const DataT* d_A_block = d_A + vbmat_idx;           //access the block on d_A.

            int k = rows_in_block, m = B_rows, n = vbmatA.block_col_size;
            int lda = rows_in_block, ldb = B_rows, ldc = C_rows;

            //multiply the blocks, store result in d_C_block
            checkCudaErrors(
                cublasGemmEx(
                    handle, CUBLAS_OP_N, CUBLAS_OP_N,
                    m,n,k,                                               //m, n, k <-- block_B: m*k   block_A: k*n   block_C: m*n
                    &alpha,
                    d_B_block,                                          // blockA device pointer,
                    data_type_AB,                                       // blockA datatype
                    ldb,                                      // blockA leading dimension
                    d_A_block,                                          // blockB device pointer
                    data_type_AB,                                       // blockB datatype
                    lda,                                             // B leading dimension
                    &beta,
                    d_C_block, data_type_C,                             // blockC device pointer, blockC type
                    ldc,                                             // C leading dimension
                    compute_type,                                       // compute_type
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

    checkCudaErrors(cublasGetMatrix(
            C_rows, C_cols, sizeof(DataT_C), d_C, C_rows, C, C_rows));

    checkCudaErrors(cudaFree(d_C));
    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_B));

    checkCudaErrors(cublasDestroy(handle));
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

void pico_print_DnM(const char* Cname, int Cn, int Cm, DataT_C* C) {

    printf("Dense matrix %s:\n", Cname);
    for (int i=0; i<Cn; i++) {
        for (int j=0; j<Cm; j++) {
            printf(" %6.3f ", C[i * Cm + j]);
        }
        printf("\n");
    }
    printf("\n");
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
        for (int i=0; i<rows; i++) {
            printf("\t\t");
            for (int j=0; j< ellValue_cols; j++) {
                if ((j%ell_blocksize) == 0)
                    printf("| ");
                if (ellColumnsInd[(i/ell_blocksize)*ellColumnsInd_cols + (j/ell_blocksize)] != -1) {
                    printf("%f ", ellValues[i*ellValue_cols +  j]);
                } else {
                    printf("Pad ");
                }
            }
            printf("\n");
            if ((i%ell_blocksize) == (ell_blocksize-1))
                printf("\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        }

        printf("\n\traw ellColInd matrix:\n\t\t");
        for (int i=0; i<ellColumnsInd_cols*ellColumnsInd_rows; i++)
            printf("%ld ", ellColumnsInd[i]);
        printf("\n");

        printf("\traw ellValue matrix:\n\t\t");
        for (int i=0; i<num_blocks*ell_blocksize*ell_blocksize; i++)
            printf("%6.3f ", ellValues[i]);
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
    int *d_RowPtr, *d_ColInd;
    DataT *d_Val;

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
        if (cmat.pattern_only == 0)
            std::copy(cmat.ma[i], cmat.ma[i] + cmat.nzcount[i], (*csrVal) + nnz);
        else
            for (int j=0; j<cmat.nzcount[i]; j++)
                (*csrVal)[nnz +j] = 1;
        nnz += cmat.nzcount[i];

    }

    return 0;
}

void cusparse_blockmat_multiplyAB(CSR& A, DataT* B, int B_cols, DataT_C* C, int C_cols, float& dt) {

    DataT *csrVal;
    int *csrRowPtr, *csrColInd;

    prepare_cusparse_CSR( A, &csrRowPtr, &csrColInd, &csrVal);
    cudaDeviceSynchronize();

    cusparse_gemm_custom(A.rows, A.cols, (int) A.nztot(), csrRowPtr, csrColInd, csrVal, B, B_cols, B_cols, C, C_cols, 1, 1, dt);

    free(csrVal);
    free(csrColInd);
    free(csrRowPtr);

    return;
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
    *ell_blocksize = A->block_col_size;
    if ((A->rows)%(*ell_blocksize) != 0 || (A->cols)%(*ell_blocksize)) {
        if ((A->rows)%(*ell_blocksize) != 0)
            printf("The number of rows is not multiple of ell_blocksize\n");
        else
            printf("The number of cols is not multiple of ell_blocksize\n");
        exit(__LINE__);
    }

    *ellColInd_rows = ((A->rows)/(*ell_blocksize));
    *ellColInd_cols = 0;
    for (int i=0; i<(*ellColInd_rows); i++) {
        if (A->nzcount[i] > (*ellColInd_cols))
            (*ellColInd_cols) = A->nzcount[i];
    }
    *ellValue_cols = (*ellColInd_cols) * (*ell_blocksize);

    *num_blocks = (*ellColInd_rows)*(*ellColInd_cols);
    *ellColInd = (intT*)malloc(sizeof(intT)*(*ellColInd_cols)*(*ellColInd_rows));
    *ellValues = (DataT_C*)malloc(sizeof(DataT)*((*num_blocks))*(*ell_blocksize)*(*ell_blocksize));

    int k_col = 0, k, i, j, vbr_blk_row_shift, bel_blk_row_shift;

    // to complete ellColInd
    for (i=0; i<(*ellColInd_rows); i++)
        for (j=0; j<(*ellColInd_cols); j++)
            if (j < A->nzcount[i]) {
                (*ellColInd)[i*(*ellColInd_cols)+j] = A->jab[k_col];
                k_col++;
            } else {
                (*ellColInd)[i*(*ellColInd_cols)+j] = -1; // '-1' means thet the algorithm autmaticly pad it with a compleate zero block
            }

    // to complete ellValues
    vbr_blk_row_shift = 0;
    bel_blk_row_shift = 0;
    for(k=0; k<(*ellColInd_rows); k++) {
        for (i=0; i<(*ell_blocksize); i++)
            for (j=0; j<(*ellValue_cols); j++)
                (*ellValues)[bel_blk_row_shift + i*(*ellValue_cols) + j] =
                    ((*ellColInd)[k*(*ellColInd_cols) + (j/(*ell_blocksize))] != -1) ? A->mab[vbr_blk_row_shift + j*(*ell_blocksize) + i] : 0.0;

        vbr_blk_row_shift += (A->nzcount[k])*(*ell_blocksize)*(*ell_blocksize);
        bel_blk_row_shift += (*ellColInd_cols)*(*ell_blocksize)*(*ell_blocksize);
    }

    return 0;
}

void bellpack_blockmat_multiplyAB(VBR* A, DataT* B, int B_cols, DataT_C* C, int C_cols, float& dt, int verbose) {

        // ellValue_cols, int *ell_blocksize, int *ellColInd_rows, int *ellColInd_cols, int *num_blocks, intT** ellColInd, DataT_C** ellValues
        int ell_blocksize, ellColInd_rows, ellColInd_cols, ellValue_cols, num_blocks;
        intT* ellColInd;
        DataT_C* ellValues;
        prepare_cusparse_BLOCKEDELLPACK(A, &ell_blocksize, &ellValue_cols, &ellColInd_rows, &ellColInd_cols, &num_blocks, &ellColInd, &ellValues);

        if (verbose > 0) {
            int pad_num = 0;
            for (int i=0; i<ellColInd_rows*ellColInd_cols; i++)
                if (ellColInd[i] == -1)
                    pad_num++;
            printf("ell_blocksize = %d, ellColInd has dimensions %d x %d with %d padding blocks, ellValues has dimensions %ld x %d\n", ell_blocksize, ellColInd_rows, ellColInd_cols, pad_num, A->rows, ellValue_cols);
        }

        cusparse_gemm_custom_ellpack(A->rows, A->cols, ell_blocksize, ellValue_cols, ellColInd_cols, ellColInd_rows, num_blocks, ellColInd, ellValues, B, B_cols, B_cols, C, C_cols, 1, 1, dt);

        free(ellColInd);
        free(ellValues);

    return;
}
