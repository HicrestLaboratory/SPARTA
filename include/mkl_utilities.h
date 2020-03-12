#include "mkl.h"

void mkl_batch_custom(const CBLAS_LAYOUT Layout,
                      const CBLAS_TRANSPOSE* transa_array, const CBLAS_TRANSPOSE* transb_array,
                      const MKL_INT* m_array, const MKL_INT* n_array, const MKL_INT* k_array,
                      const float* alpha_array,
                      const float **a_array, const MKL_INT* lda_array,
                      const float **b_array, const MKL_INT* ldb_array,
                      const float* beta_array,
                      float **c_array, const MKL_INT* ldc_array,
                      const MKL_INT group_count, const MKL_INT* group_size);

void mkl_batch_custom(const CBLAS_LAYOUT Layout,
                      const CBLAS_TRANSPOSE* transa_array, const CBLAS_TRANSPOSE* transb_array,
                      const MKL_INT* m_array, const MKL_INT* n_array, const MKL_INT* k_array,
                      const double* alpha_array,
                      const double **a_array, const MKL_INT* lda_array,
                      const double **b_array, const MKL_INT* ldb_array,
                      const double* beta_array,
                      double **c_array, const MKL_INT* ldc_array,
                      const MKL_INT group_count, const MKL_INT* group_size);

void mkl_gemm_custom(const CBLAS_LAYOUT Layout,
                  const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb,
                  const MKL_INT m, const MKL_INT n, const MKL_INT k,
                  const float alpha, const float *a, const MKL_INT lda,
                  const float *b, const MKL_INT ldb,
                     const float beta, float *c, const MKL_INT ldc);

void mkl_gemm_custom(const CBLAS_LAYOUT Layout,
                     const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb,
                     const MKL_INT m, const MKL_INT n, const MKL_INT k,
                     const double alpha, const double *a, const MKL_INT lda,
                     const double *b, const MKL_INT ldb,
                     const double beta, double *c, const MKL_INT ldc);

void mkl_spmm_custom(const sparse_operation_t operation,
                     const float alpha, const sparse_matrix_t A,
                     const struct matrix_descr descr,
                     const sparse_layout_t layout,
                     const float *x, const MKL_INT columns,
                     const MKL_INT ldx,
                     const float beta,
                     float *y, const MKL_INT ldy);

void mkl_spmm_custom(const sparse_operation_t operation,
                     const double alpha, const sparse_matrix_t A,
                     const struct matrix_descr descr,
                     const sparse_layout_t layout,
                     const double *x, const MKL_INT columns,
                     const MKL_INT ldx,
                     const double beta,
                     double *y, const MKL_INT ldy);

void mkl_create_csr_custom(sparse_matrix_t *A,
                       const sparse_index_base_t indexing,
                       const MKL_INT rows, const MKL_INT cols,
                       MKL_INT *rows_start, MKL_INT *rows_end,
                       MKL_INT *col_indx,
                           float *values);

void mkl_create_csr_custom(sparse_matrix_t *A,
                           const sparse_index_base_t indexing,
                           const MKL_INT rows, const MKL_INT cols,
                           MKL_INT *rows_start, MKL_INT *rows_end,
                           MKL_INT *col_indx,
                           double *values);

void convert_to_MKL(SparMat &spmt, sparse_matrix_t &A);

void mkl_blockmat_multiply(const VBSparMat &VBMat,
                        DataT *X, int X_cols,
                        DataT *Y);

void mkl_blockmat_batch_multiply(const VBSparMat &VBMat,
                                 DataT *X,int X_cols,
                                 DataT *Y);






