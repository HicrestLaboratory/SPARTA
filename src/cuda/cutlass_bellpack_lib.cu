#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>

#include <cutlass/layout/matrix.h>
#include "cutlass/cutlass.h"
#include "cutlass/gemm/gemm.h"
#include "cutlass/gemm/kernel/gemm_grouped.h"
#include "cutlass/gemm/kernel/default_gemm_grouped.h"
#include "cutlass/gemm/device/ell_gemm.h"

#include "cutlass/util/tensor_view_io.h"
#include "cutlass/util/host_tensor.h"
#include "cutlass/util/reference/host/gemm.h"

#include "cutlass_bellpack_lib.h"

#define DataT float
#define DataT_Cutlass cutlass::half_t
#define NOELL_TEST

template<typename iT, typename T>
int compute_cutlass_bellpack (int rows, int cols, int ell_blocksize, int ellValue_cols, iT* ellColInd, T *ellValues, int B_rows, int B_cols, T *B_vals, int C_rows, int C_cols, T *C_vals) {

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

    printf("-------------------------------------------------------------------------------------------------\n");

    // Configure the GEMM arguments
    float alpha=1.0, beta=1.0;

    DataT_Cutlass *ptrD = tensorC.device_data();

    int lda = tensor.device_ref().stride(0);
    int ldb = tensorB.device_ref().stride(0);
    int ldc = tensorC.device_ref().stride(0);
    int ldd = tensorC.device_ref().stride(0);

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

    printf("File %s at line %d\n", __FILE__, __LINE__);

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
