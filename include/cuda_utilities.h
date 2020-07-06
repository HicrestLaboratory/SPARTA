#include "comp_mats.h"

void cublas_blockmat_multiply(const VBS& vbmatA, float* B, int B_cols, int B_lead_dim, float* C, int C_lead_dim);

int cublas_gemm_custom(const float* A, unsigned int A_rows, unsigned int A_cols, unsigned int lda, const float* B, unsigned int B_cols, unsigned int ldb, float* C, unsigned int ldc, const float alpha, const float beta);

int cusparse_gemm_custom(int rows, int cols, int nnz, int* csrRowPtr, int* csrColInd, float* csrVal, float* B, int B_cols, int B_lead_dim, float* C, int C_lead_dim, const float alpha, const float beta);

int prepare_cusparse_CSR(CSR& cmat, int **csrRowPtr, int **csrColInd, float **csrVal);


