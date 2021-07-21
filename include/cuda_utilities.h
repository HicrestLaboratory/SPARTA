#include "comp_mats.h"

void cublas_blockmat_multiply(const VBS& vbmatA, DataT* B, int B_cols, int B_lead_dim, DataT* C, int C_lead_dim, float &dt, int n_streams);

int cublas_gemm_custom(const DataT* A, unsigned int A_rows, unsigned int A_cols, unsigned int lda, const DataT* B, unsigned int B_cols, unsigned int ldb, DataT* C, unsigned int ldc, const DataT alpha, const DataT beta, float& dt);

int cusparse_gemm_custom(int rows, int cols, int nnz, int* csrRowPtr, int* csrColInd, float* csrVal, float* B, int B_cols, int B_lead_dim, float* C, int C_lead_dim, const float alpha, const float beta, float& dt);

int prepare_cusparse_CSR(CSR& cmat, int *csrRowPtr, int *csrColInd, float *csrVal);


