// Utilities and system includes
#include <assert.h>
#include <helper_string.h>  // helper for shared functions common to CUDA Samples

// CUDA runtime
#include <cuda_runtime.h>
#include <cublas_v2.h>

// CUDA and CUBLAS functions
#include <helper_functions.h>
#include <helper_cuda.h>

//Matrix-Matrix multiplication with cublas. A,B,C are in column-major order.
int cublas_gemm_custom(float *A, unsigned int uiWA, unsigned int uiHA, unsigned int lda,
                   float *B, unsigned int uiWB, unsigned int ldb,
                   float *C, unsigned int ldc)
{
    int block_size = 32;
    
    //deduce matrices dimensions
    unsigned int uiHB = uiWA;
    unsigned int uiWC = uiWB;
    unsigned int uiHC = uiHA;
    
    unsigned int size_A = uiWA * uiHA;
    unsigned int mem_size_A = sizeof(float) * size_A;

    unsigned int size_B = uiWB * uiHB;
    unsigned int mem_size_B = sizeof(float) * size_B;
    
    unsigned int size_C = uiWC * uiHC;
    unsigned int mem_size_C = sizeof(float) * size_C;
    
    // allocate device memory
    float *d_A, *d_B, *d_C;
    checkCudaErrors(cudaMalloc((void **) &d_A, mem_size_A));
    checkCudaErrors(cudaMalloc((void **) &d_B, mem_size_B));
    checkCudaErrors(cudaMalloc((void **) &d_C, mem_size_C));
    
    //copy matrices to device
    checkCudaErrors(cublasSetMatrix(
                                    uiHA,uiWA, sizeof(float), A, lda, d_A, int ldb));
    checkCudaErrors(cublasSetMatrix(
                                    uiHB,uiWB, sizeof(float), B, ldb, d_B, int ldb));

    // setup execution parameters
    dim3 threads(block_size, block_size);
    dim3 grid(uiWC / threads.x, uiHC / threads.y);

    // CUBLAS version 2.0
    const float alpha = 1.0f;
    const float beta  = 1.0f;
    cublasHandle_t handle;

    checkCudaErrors(cublasCreate(&handle));

    //Perform gemm operation with cublas
    checkCudaErrors(cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                                uiWB, uiHA, uiWA,
                                &alpha,
                                d_A, lda,
                                d_B, ldb,
                                &beta,
                                d_C, uiHC));

    // copy result from device to host
    checkCudaErrors(cublasGetMatrix(uiHC,uiWC, sizeof(float), d_C, uiHC, C, ldc));

    // Destroy the handle
    checkCudaErrors(cublasDestroy(handle));

    // clean up memory
    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_B));
    checkCudaErrors(cudaFree(d_C));

    if (resCUBLAS == true)
    {
        return EXIT_SUCCESS;    // return value = 1
    }
    else
    {
        return EXIT_FAILURE;     // return value = 0
    }
}


void randomInit(float *data, int size)
{
    for (int i = 0; i < size; ++i)
        data[i] = rand() / (float)RAND_MAX;
}

int main(int argc, char **argv){
    int AW = 10;
    int AH = 5;
    int sizeA = AW*AH;
    float A[sizeA];
    
    int BW = 8;
    int BH = AW;
    int sizeB = BW*BH;
    float B[sizeB];
    
    int CH = AH;
    int CW = BW;
    int sizeC = CW*CH;
    float C[sizeC];
    
    randomInit(A,sizeA);
    randomInit(B,sizeB);
    
    return cublas_gemm_custom(A, AW, AH, AH,
                       B, BW, BH,
                       C, CH);
}
