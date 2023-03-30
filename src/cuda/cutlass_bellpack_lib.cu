#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>

//CUDA Utilities and system includes
#include <assert.h>

// CUDA runtime
#include <cuda_runtime.h>
#include "helper_cuda.h"

#include <cutlass/layout/matrix.h>
#include "cutlass/cutlass.h"
#include "cutlass/gemm/gemm.h"
#include "cutlass/gemm/kernel/gemm_grouped.h"
#include "cutlass/gemm/kernel/default_gemm_grouped.h"
#include "cutlass/gemm/device/ell_gemm.h"

#include "cutlass/util/tensor_view_io.h"
// #include "cutlass/util/host_tensor.h"
#include "cutlass/util/reference/host/gemm.h"

// --------------
#include "cutlass/gemm/device/gemm_sparse.h"
#include "cutlass/util/host_reorder.h"
#include "cutlass/util/host_uncompress.h"
#include "cutlass/util/reference/host/tensor_compare.h"
#include "cutlass/util/reference/host/tensor_copy.h"
#include "helper.h"

#include <cutlass/numeric_types.h>
#include <cutlass/gemm/device/gemm.h>

#include <cutlass/util/host_tensor.h>
// --------------

#include "cuda_utilities.h"
#include "cutlass_bellpack_lib.h"


#define DataT float
#define DataT_Cutlass cutlass::half_t
#define NOELL_TEST

template<typename iT, typename T>
int compute_cutlass_bellpack (int rows, int cols, int ell_blocksize, int ellValue_cols, iT* ellColInd, T *ellValues, int B_rows, int B_cols, T *B_vals, int C_rows, int C_cols, T *C_vals, float& dt) {

    int i, j, ellColInd_rows = (rows/ell_blocksize), ellColInd_cols = (ellValue_cols/ell_blocksize);

//     printf("Inside of %s\n", __func__);
//     print_bellpack("Bellpack_A", rows, cols, ell_blocksize, ellValue_cols, ellColInd_rows, ellColInd_cols, ellColInd_rows*ellColInd_cols, ellColInd, ellValues);

    cutlass::HostTensor<DataT_Cutlass, cutlass::layout::RowMajor> tensor({ellColInd_rows*ell_blocksize, ellColInd_cols*ell_blocksize});
    for (i = 0; i < (ellColInd_rows*ell_blocksize); ++i) {
        for (j = 0; j < (ellColInd_cols*ell_blocksize); ++j) {
            // Write the element at location {i, j} in host memory
            tensor.host_ref().at({i, j}) = (DataT_Cutlass)ellValues[i*ellValue_cols + j];
        }
    }
    // Copy host memory to device memory
    tensor.sync_device();
    // Obtain a device pointer usable in CUDA kernels
    DataT_Cutlass *device_ptr = tensor.device_data();


    cutlass::HostTensor<int32_t, cutlass::layout::RowMajor> tensor_ellIdx({ellColInd_rows, ellColInd_cols});
    for (i = 0; i < ellColInd_rows; ++i) {
        for (j = 0; j < ellColInd_cols; ++j) {
            // Write the element at location {i, j} in host memory
            tensor_ellIdx.host_ref().at({i, j}) = (int32_t)ellColInd[i*ellColInd_cols + j];
        }
    }
    // Copy host memory to device memory
    tensor_ellIdx.sync_device();
    // Obtain a device pointer usable in CUDA kernels
    int32_t *device_ptr_ellIdx = tensor_ellIdx.device_data();

//     printf("Inside of %s\n", __func__);
//     std::cout << "view of tensor_ellIdx:" << std::endl;
//     cutlass::TensorView<int32_t, cutlass::layout::RowMajor> view_ellIdx = tensor_ellIdx.host_view();
//     std::cout << view_ellIdx << std::endl;
//     std::cout << std::endl;


    cutlass::HostTensor<DataT_Cutlass, cutlass::layout::ColumnMajor> tensorB({B_rows, B_cols});
    for (i = 0; i < (B_rows); ++i) {
        for (j = 0; j < (B_cols); ++j) {
            // Write the element at location {i, j} in host memory
            tensorB.host_ref().at({i, j}) = (DataT_Cutlass)B_vals[i*B_cols + j];
        }
    }
    // Copy host memory to device memory
    tensorB.sync_device();
    // Obtain a device pointer usable in CUDA kernels
    DataT_Cutlass *device_ptrB = tensorB.device_data();

//     printf("Inside of %s\n", __func__);
//     std::cout << "view of tensorB:" << std::endl;
//     cutlass::TensorView<DataT_C, cutlass::layout::ColumnMajor> viewB = tensorB.host_view();
//     std::cout << viewB << std::endl;
//     std::cout << std::endl;


    cutlass::HostTensor<DataT_Cutlass, cutlass::layout::ColumnMajor> tensorC({rows, B_cols});
    for (i = 0; i < (rows); ++i) {
        for (j = 0; j < (B_cols); ++j) {
            // Write the element at location {i, j} in host memory
            tensorC.host_ref().at({i, j}) = (DataT_Cutlass) 0.0;
        }
    }
    // Copy host memory to device memory
    tensorC.sync_device();
    // Obtain a device pointer usable in CUDA kernels
    DataT_Cutlass *device_ptrC = tensorC.device_data();

    // ================================================================================

    cudaDeviceProp props;

    cudaError_t error = cudaGetDeviceProperties(&props, 0);
    if (error != cudaSuccess) {
        std::cerr << "cudaGetDeviceProperties() returned an error: " << cudaGetErrorString(error) << std::endl;
        return -1;
    }

    if (__CUDACC_VER_MAJOR__ < 11 || props.major < 8) {

        //
        // This example requires an NVIDIA Ampere-architecture GPU.
        //

        std::cout
        << "CUTLASS's BlockedEll SpMM example requires a GPU of NVIDIA's Ampere Architecture or "
        << "later (compute capability 80 or greater).\n";
    }


    //
    // Define the BlockedEll type
    //

    using Gemm = typename cutlass::gemm::device::EllGemm<
        DataT_Cutlass,
        cutlass::layout::RowMajor,
        DataT_Cutlass,
        cutlass::layout::ColumnMajor,
        DataT_Cutlass,
        cutlass::layout::ColumnMajor,
        float,
        cutlass::arch::OpClassTensorOp,
        cutlass::arch::Sm80>;
    Gemm gemm_op;
    cutlass::Status status;

// -------------------------------------------------------------------------------------------------

    // Configure the GEMM arguments
    float alpha=1.0, beta=1.0;

    DataT_Cutlass *ptrD = tensorC.device_data();

    int lda = tensor.device_ref().stride(0);
    int ldb = tensorB.device_ref().stride(0);
    int ldc = tensorC.device_ref().stride(0);
    int ldd = tensorC.device_ref().stride(0);

    //initialize cuda events
    cudaEvent_t start, stop;
    checkCudaErrors(cudaEventCreate(&start));
    checkCudaErrors(cudaEventCreate(&stop));
    checkCudaErrors(cudaEventRecord(start, 0));

    // Configure GEMM arguments
    status = gemm_op({
      {rows, B_cols, cols},
      {tensor.device_ref(), lda},
      {tensorB.device_ref(), ldb},
      {tensorC.device_ref(), ldc},
      {ptrD, ldd},
      tensor_ellIdx.device_data(),
      ellValue_cols,
      ell_blocksize,
      0 /*options.a_base*/,
      {alpha, beta}
    });

    if (status != cutlass::Status::kSuccess) {
      std::cerr << "Failed to initialize CUTLASS BlockedEll SpMM kernel." << std::endl;
      return(__LINE__);
    }

    if (status != cutlass::Status::kSuccess) {
      std::cerr << "Failed to run CUTLASS BlockedEll SpMM kernel." << std::endl;
      return(__LINE__);
    }

    // Wait for completion
    error = cudaDeviceSynchronize();
    checkCudaErrors( cudaEventRecord(stop, 0) );
    checkCudaErrors( cudaEventSynchronize(stop) );
    checkCudaErrors( cudaEventElapsedTime(&dt, start, stop) );
    checkCudaErrors( cudaEventDestroy(start) );
    checkCudaErrors( cudaEventDestroy(stop) );

    if (error != cudaSuccess)  {
      std::cerr << "Kernel execution error: " << cudaGetErrorString(error);
      return(__LINE__);
    }

    // ================================================================================

    tensorC.sync_host();
//     printf("Inside of %s\n", __func__);
//     std::cout << "view of tensorC:" << std::endl;
//     cutlass::TensorView<DataT_C, cutlass::layout::ColumnMajor> viewC = tensorC.host_view();
//     std::cout << viewC << std::endl;
//     std::cout << std::endl;
    for (i = 0; i < C_rows; ++i) {
        for (j = 0; j < C_cols; ++j) {
            // Write the element at location {i, j} in host memory
            C_vals[i*C_cols + j] = (T)tensorC.host_ref().at({i, j});
        }
    }

    return(0);

}

void bellpack_cutlass_multiplyAB(VBR* A, DataT* B, int B_cols, DataT_C* C, int C_cols, float& dt, int verbose) {

        // ellValue_cols, int *ell_blocksize, int *ellColInd_rows, int *ellColInd_cols, int *num_blocks, intT** ellColInd, Cutlass_DataT** ellValues
        int ell_blocksize, ellColInd_rows, ellColInd_cols, ellValue_cols, num_blocks;
        intT* ellColInd;
        DataT* ellValues;
        prepare_cusparse_BLOCKEDELLPACK(A, &ell_blocksize, &ellValue_cols, &ellColInd_rows, &ellColInd_cols, &num_blocks, &ellColInd, &ellValues);

        if (verbose > 1) {
            int pad_num = 0;
            for (int i=0; i<ellColInd_rows*ellColInd_cols; i++)
                if (ellColInd[i] == -1)
                    pad_num++;
            printf("ell_blocksize = %d, ellColInd has dimensions %d x %d with %d padding blocks, ellValues has dimensions %ld x %d\n", ell_blocksize, ellColInd_rows, ellColInd_cols, pad_num, A->rows, ellValue_cols);
        }

        compute_cutlass_bellpack<intT,DataT>(A->rows, A->cols, ell_blocksize, ellValue_cols, ellColInd, ellValues, A->cols, B_cols, B, A->rows, B_cols, C, dt);

        free(ellColInd);
        free(ellValues);

    return;
}

int cutlass_dense_multiplyAB(int m, int k, DataT* inputA, int n, DataT* inputB, float alp, float bet, DataT_C* output, float& dt) {

  dt = -1.0;

  // Define the GEMM operation
  using Gemm = cutlass::gemm::device::Gemm<
    cutlass::half_t,                           // ElementA
    cutlass::layout::ColumnMajor,              // LayoutA
    cutlass::half_t,                           // ElementB
    cutlass::layout::ColumnMajor,              // LayoutB
    cutlass::half_t,                           // ElementOutput
    cutlass::layout::ColumnMajor,              // LayoutOutput
    float,                                     // ElementAccumulator
    cutlass::arch::OpClassTensorOp,            // tag indicating Tensor Cores
    cutlass::arch::Sm80                        // tag indicating target GPU compute architecture
  >;

  Gemm gemm_op;
  cutlass::Status status;

  //
  // Define the problem size
  //
  int M = m;
  int N = n;
  int K = k;

  float alpha = alp;
  float beta  = bet;

  //
  // Allocate device memory
  //



  cutlass::HostTensor<cutlass::half_t, cutlass::layout::ColumnMajor> A({M, K});
  cutlass::HostTensor<cutlass::half_t, cutlass::layout::ColumnMajor> B({K, N});
  cutlass::HostTensor<cutlass::half_t, cutlass::layout::ColumnMajor> C({M, N});

  // Input matrices to cutlass' structures



  for (int i=0; i<M; ++i)
    for (int j=0; j<K; ++j)
      A.host_ref().at({i, j}) = inputA[i*K +j];
  // Copy host memory to device memory
  A.sync_device();



  for (int i=0; i<K; ++i)
    for (int j=0; j<N; ++j)
      B.host_ref().at({i, j}) = inputB[i*N +j];
  // Copy host memory to device memory
  B.sync_device();



  cutlass::half_t const *ptrA = A.device_data();
  cutlass::half_t const *ptrB = B.device_data();
  cutlass::half_t const *ptrC = C.device_data();
  cutlass::half_t       *ptrD = C.device_data();



  int lda = A.device_ref().stride(0);
  int ldb = B.device_ref().stride(0);
  int ldc = C.device_ref().stride(0);
  int ldd = C.device_ref().stride(0);

  //
  // Launch GEMM on the device
  //



  status = gemm_op({
    {M, N, K},
    {ptrA, lda},            // TensorRef to A device tensor
    {ptrB, ldb},            // TensorRef to B device tensor
    {ptrC, ldc},            // TensorRef to C device tensor
    {ptrD, ldd},            // TensorRef to D device tensor - may be the same as C
    {alpha, beta}           // epilogue operation arguments
  });



  if (status != cutlass::Status::kSuccess) {
    return -1;
  }



  // Copy host memory to device memory
  C.sync_host();
  for (int i=0; i<M; ++i)
    for (int j=0; j<N; ++j)
      output[i*N +j] = C.host_ref().at({i, j});



  return 0;
}

// void cutlass_dense_multiplyAB(int rows, int cols, DataT* A, int B_rows, int B_cols, DataT* B, DataT_C* C, float& dt) {
//     int i, j;
//     cutlass::HostTensor<DataT_Cutlass, cutlass::layout::ColumnMajor> tensor({rows, cols});
//     for (i = 0; i < rows; ++i) {
//         for (j = 0; j < cols; ++j) {
//
//             // Write the element at location {i, j} in host memory
//             tensor.host_ref().at({i, j}) = (DataT_Cutlass)A[i*cols + j];
//
//         }
//     }
//
//     // Copy host memory to device memory
//     tensor.sync_device();
//
//     // Obtain a device pointer usable in CUDA kernels
//     DataT_Cutlass *device_ptr = tensor.device_data();
//
// //     std::cout << "view of tensor:" << std::endl;
// //     cutlass::TensorView<DataT_Cutlass, cutlass::layout::ColumnMajor> view = tensor.host_view();
// //     std::cout << view << std::endl;
// //     std::cout << std::endl;
//
//
//     cutlass::HostTensor<DataT_Cutlass, cutlass::layout::ColumnMajor> tensorB({B_rows, B_cols});
//
//
//     for (i = 0; i < (B_rows); ++i) {
//         for (j = 0; j < (B_cols); ++j) {
//
//             // Write the element at location {i, j} in host memory
//             tensorB.host_ref().at({i, j}) = B[i*B_cols + j];
//
//         }
//     }
//
//     // Copy host memory to device memory
//     tensorB.sync_device();
//
//     // Obtain a device pointer usable in CUDA kernels
//     DataT_Cutlass *device_ptrB = tensorB.device_data();
//
// //     std::cout << "view of tensorB:" << std::endl;
// //     cutlass::TensorView<DataT_Cutlass, cutlass::layout::ColumnMajor> viewB = tensorB.host_view();
// //     std::cout << viewB << std::endl;
// //     std::cout << std::endl;
//
//     cutlass::HostTensor<DataT_Cutlass, cutlass::layout::ColumnMajor> tensorC({rows, B_cols});
//
//
//     for (i = 0; i < (rows); ++i) {
//         for (j = 0; j < (B_cols); ++j) {
//
//             // Write the element at location {i, j} in host memory
//             tensorC.host_ref().at({i, j}) = (DataT_Cutlass) 0.0;
//
//         }
//     }
//
//     // Copy host memory to device memory
//     tensorC.sync_device();
//
//     // Obtain a device pointer usable in CUDA kernels
//     DataT_Cutlass *device_ptrC = tensorC.device_data();
//
// //     std::cout << "view of tensorC:" << std::endl;
// //     cutlass::TensorView<DataT_Cutlass, cutlass::layout::ColumnMajor> viewC = tensorC.host_view();
// //     std::cout << viewC << std::endl;
// //     std::cout << std::endl;
//
//   // ------------------------------------------------------------------------------
//
//   using Gemm = cutlass::gemm::device::Gemm<
//     DataT_Cutlass,                           // ElementA
//     cutlass::layout::ColumnMajor,              // LayoutA
//     DataT_Cutlass,                           // ElementB
//     cutlass::layout::ColumnMajor,              // LayoutB
//     DataT_Cutlass,                           // ElementOutput
//     cutlass::layout::ColumnMajor,              // LayoutOutput
//     float,                                     // ElementAccumulator
//     cutlass::arch::OpClassTensorOp,            // tag indicating Tensor Cores
//     cutlass::arch::Sm80                        // tag indicating target GPU compute architecture
//   >;
//
//   Gemm gemm_op_dn;
//   cutlass::Status status;
//
//   //
//   // Launch GEMM on the device
//   //
//
//   float alpha=1.0, beta=1.0;
//   int m = rows, k = cols, n = B_cols;
//
//   DataT_Cutlass *ptrD = tensorC.device_data();
//
//   int lda = tensor.device_ref().stride(0);
//   int ldb = tensorB.device_ref().stride(0);
//   int ldc = tensorC.device_ref().stride(0);
//   int ldd = tensorC.device_ref().stride(0);
//
//   //initialize cuda events
//   cudaEvent_t start, stop;
//   checkCudaErrors(cudaEventCreate(&start));
//   checkCudaErrors(cudaEventCreate(&stop));
//   checkCudaErrors(cudaEventRecord(start, 0));
//
//   status = gemm_op_dn({
//     {m, n, k},
//     {device_ptr, lda},
//     {device_ptrB, ldb},
//     {device_ptrC, ldc},
//     {ptrD, ldd},
//     {alpha, beta}
//   });
//
//   if (status != cutlass::Status::kSuccess) {
//     fprintf(stderr, "ERROR in file %s at line %d\n", __FILE__, __LINE__);
//   } else {
//     printf("SUCCESS in file %s at line %d\n", __FILE__, __LINE__);
//   }
//
//   // Wait for completion
//   checkCudaErrors( cudaDeviceSynchronize() );
//   checkCudaErrors( cudaEventRecord(stop, 0) );
//   checkCudaErrors( cudaEventSynchronize(stop) );
//   checkCudaErrors( cudaEventElapsedTime(&dt, start, stop) );
//   checkCudaErrors( cudaEventDestroy(start) );
//   checkCudaErrors( cudaEventDestroy(stop) );
//
//   tensorC.sync_host();
// //   std::cout << "view of tensorC:" << std::endl;
// //   std::cout << viewC << std::endl;
// //   std::cout << std::endl;
//
//   for (i = 0; i < (rows); ++i) {
//         for (j = 0; j < (B_cols); ++j) {
//
//             // Write the element at location {i, j} in host memory
//             C[i * B_cols + j] = (DataT_C)tensorC.host_ref().at({i, j});
//
//         }
//     }
//
//   return;
// }
