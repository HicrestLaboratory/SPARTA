#include "matrices.h"

void cublas_blockmat_multiply(const VBR& vbmatA, DataT* B, int B_rows, int B_lead_dim, DataT_C* C, int C_lead_dim, float &dt, int n_streams = 16);

//int cublas_gemm_custom(const DataT* A, unsigned int A_rows, unsigned int A_cols, unsigned int lda, const DataT* B, unsigned int B_cols, unsigned int ldb, DataT_C* C, unsigned int ldc, const DataT_C alpha, const DataT_C beta, float& dt);

//int cusparse_gemm_custom(int rows, int cols, int nnz, int* csrRowPtr, int* csrColInd, DataT* csrVal, DataT* B, int B_cols, int B_lead_dim, DataT_C* C, int C_lead_dim, const DataT_C alpha, const DataT_C beta, float& dt);

//int prepare_cusparse_CSR(CSR& cmat, int *csrRowPtr, int *csrColInd, DataT*csrVal);


