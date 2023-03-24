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

#define intT int32_t
#define DataT float
#define DataT_Cutlass cutlass::half_t
#define NOELL_TEST

void print_bellpack(const char* Aname, int rows, int cols, int ell_blocksize, int ellValue_cols, int ellColumnsInd_rows, int ellColumnsInd_cols, int num_blocks, intT* ellColumnsInd, DataT* ellValues) {

    printf("Blocked-ellpac matrix %s:\n", Aname);
    printf("\tellColInd matrix:\n");
    for (int i=0; i<ellColumnsInd_rows; i++) {
        printf("\t\t");
        for (int j=0; j< ellColumnsInd_cols; j++)
            printf("%d ", ellColumnsInd[i*ellColumnsInd_cols + j]);
        printf("\n");
    }

    printf("\tellValue matrix:\n");
    if (rows < 50 && ellValue_cols < 50) {
        for (int i=0; i<rows; i++) {
            printf("\t\t");
            for (int j=0; j< ellValue_cols; j++) {
                if ((j%ell_blocksize) == 0)
                    printf("| ");
                if (ellColumnsInd[(i/ell_blocksize)*ellColumnsInd_cols + (j/ell_blocksize)] != -1) {
                    printf("%4.3f ", (double) ellValues[i*ellValue_cols +  j]);
                } else {
                    printf("Pad ");
                }
            }
            printf("\n");
            if ((i%ell_blocksize) == (ell_blocksize-1))
                printf("\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        }

        printf("\traw ellValue matrix:\n\t\t");
        for (int i=0; i<num_blocks*ell_blocksize*ell_blocksize; i++)
            printf("%6.3f ", (double) ellValues[i]);
        printf("\n\n");
    } else {
        for (int i=0; i<ellColumnsInd_rows; i++) {
            printf("\t\t");
            for (int j=0; j< ellColumnsInd_cols; j++) {
                if (ellColumnsInd[i*ellColumnsInd_cols + j] != -1) {
                    printf("|%5.1f-", (double)ellValues[(i*ell_blocksize)*ellValue_cols + (j*ell_blocksize)]);
                    printf("%5.1f", (double)ellValues[((i+1)*ell_blocksize-1)*ellValue_cols + ((j+1)*ell_blocksize-1)]);
                } else {
                    printf("| Padding |");
                }
            }
            printf("|\n");
        }
    }

    if (rows < 50 && ellValue_cols < 50) {
        printf("\n\traw ellColInd matrix:\n\t\t");
        for (int i=0; i<ellColumnsInd_cols*ellColumnsInd_rows; i++)
            printf("%d ", ellColumnsInd[i]);
        printf("\n");
    }

}

void print_dense(const char* Aname, int rows, int cols, DataT* values) {

    printf("Dense matrix %s:\n", Aname);
    for (int i=0; i<rows; i++) {
        for (int j=0; j<cols; j++)
            printf("%6.3f ", (double) values[i*cols + j]);
        printf("\n");
    }
    printf("\n\n");
}

void print_dense(const char* Aname, int rows, int cols, int* values) {

    printf("Dense matrix %s:\n", Aname);
    for (int i=0; i<rows; i++) {
        for (int j=0; j<cols; j++)
            printf("%3d ", values[i*cols + j]);
        printf("\n");
    }
    printf("\n\n");
}

float gaussSum (float start, float stop, float increment) {
    float len = (stop - start)/increment;
    return (((start+stop)*(stop-start+1))/((float)2));
}

int gen_bellpack_sumGaussMat (int rows, int cols, int ell_blocksize, int ellValue_cols, intT* ellColumnsInd, DataT** ellValues, DataT** resvect) {
    int i, j;
    float k, increment = 0.5;

    // Check inputs
    if ((rows % ell_blocksize != 0) && (cols % ell_blocksize != 0)) {
        printf("rows and/or cols is not multiple of ell_blocksize\n");
        return (1);
    }
    int ellColInd_rows = rows / ell_blocksize, ellColInd_cols = ellValue_cols / ell_blocksize;
    for (i=0; i<(ellColInd_rows*ellColInd_cols); i++)
        if (ellColumnsInd[i] >= (cols/ell_blocksize)) {
            printf("ellColumnsInd[%d] >= (cols/ell_blocksize)\n", i);
            return(1);
        }

    // Generate solution vector
    intT nnz_per_blkline[ellColInd_rows];
    for (i=0; i<ellColInd_rows; i++) {
        nnz_per_blkline[i] = 0;
        for (j=0; j<ellColInd_cols; j++)
            if (ellColumnsInd[i*ellColInd_cols + j] != -1)
                nnz_per_blkline[i]++;
            else
                j = ellColInd_cols;
    }

    float nnz_per_line[rows];
    for (i=0; i<rows; i++)
        nnz_per_line[i] = (float)(nnz_per_blkline[i/ell_blocksize]);

    // Generate ellValues matrix
    (*ellValues) = (DataT*) malloc(sizeof(DataT)*rows*ellValue_cols);
    k = 0.0;
    for (i=0; i<rows; i++)
        for (j=0; j<ellValue_cols; j++) {
//             printf("(%d, %d) ---> (%d, %d) ==> %d\n", i, j, i/ell_blocksize, j/ell_blocksize, ellColumnsInd[(i/ell_blocksize)*ellColInd_cols + (j/ell_blocksize)]);
            if (ellColumnsInd[(i/ell_blocksize)*ellColInd_cols + (j/ell_blocksize)] != -1) {
                (*ellValues)[i*ellValue_cols + j] = ((float)i)/(nnz_per_line[i]*ell_blocksize);// k;
                k += increment;
            } else {
                (*ellValues)[i*ellValue_cols + j] = 0.0;
            }
        }


    float start = 0.0, stop = (nnz_per_line[0]*ell_blocksize -1)*increment;
    (*resvect) = (DataT*) malloc(sizeof(DataT)*rows);
    (*resvect)[0] = (DataT) gaussSum(start, stop, increment);
    for (i=1; i<rows; i++) {
        start += (nnz_per_line[i-1]*ell_blocksize)*increment;
        stop  += ((nnz_per_line[i]*ell_blocksize)-1)*increment;
//         printf("start = %f, stop = %f\n", start, stop);
        (*resvect)[i] = (DataT) gaussSum(start, stop, increment);
    }

    return(0);
}

template<typename T>
int compute_cutlass_bellpack (int rows, int cols, int ell_blocksize, int ellValue_cols, intT* ellColInd, T *ellValues, int B_rows, int B_cols, DataT *B_vals, int C_rows, int C_cols, T *C_vals) {

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

int main(void) {

    int rows = 128, cols = 256, ell_blocksize = 16, ellValue_cols = 128, i, j, k, tmp;
    intT* ellColInd;
    DataT *ellValues, *resvect;

//     prepare_cusparse_BLOCKEDELLPACK(A, &ell_blocksize, &ellValue_cols, &ellColInd_rows, &ellColInd_cols, &num_blocks, &ellColInd, &ellValues);
    ellColInd = (intT*) malloc(sizeof(intT) * (rows/ell_blocksize) * (ellValue_cols/ell_blocksize));
    for (i=0; i<(rows/ell_blocksize); i++) {
        k=0;
        for (j=0; j<(ellValue_cols/ell_blocksize); j++) {
            tmp=0;
            while (tmp==0) {
                if (i==0 || j == 0 || (k>=(cols/ell_blocksize)) || ( rand() < (RAND_MAX/3) )) {
                    ellColInd[i*(ellValue_cols/ell_blocksize)+j] = ((k<(cols/ell_blocksize)) ? k : -1);
                    tmp = 1;
                }
                k++;
            }
        }
    }

    gen_bellpack_sumGaussMat (rows, cols, ell_blocksize, ellValue_cols, ellColInd, &ellValues, &resvect);

    int ellColInd_rows = rows/ell_blocksize, ellColInd_cols = ellValue_cols/ell_blocksize, num_blocks = ellColInd_rows*ellColInd_cols;
    print_bellpack("Bellpack_A", rows, cols, ell_blocksize, ellValue_cols, ellColInd_rows, ellColInd_cols, num_blocks, ellColInd, ellValues);

//     printf("result vector:\n\t");
//     for (i=0; i<rows; i++)
//         printf("%f ", (float) resvect[i]);
//     printf("\n\n");

    int B_rows = cols, B_cols = 3;
    DataT *B_vals = (DataT*) malloc(sizeof(DataT) * B_rows * B_cols);
    for (i=0; i<(B_rows*B_cols); i++)
        B_vals[i] = (DataT) 1.0;

    if (B_rows<50 && B_cols < 50) {
        print_dense("DnB", B_rows, B_cols, B_vals);
    }

    printf("=================================================================================================\n");

    DataT* C_vals = (DataT*) malloc(sizeof(DataT)*rows*B_cols);
    compute_cutlass_bellpack<DataT>(rows, cols, ell_blocksize, ellValue_cols, ellColInd, ellValues, B_rows, B_cols, B_vals, rows, B_cols, C_vals);
    free(ellColInd);
    free(ellValues);

    printf("compute_cutlass_bellpack DONE!\n");

    print_dense("C", rows, B_cols, C_vals);

    return(0);
}
