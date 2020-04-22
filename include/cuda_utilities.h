#include "globheads.h"
#include "protos.h"

void cublas_blockmat_multiply(const VBSparMat &VBMat, float *X, int X_cols, float *Y);

int cublas_gemm_custom(const float *A, unsigned int A_cols, unsigned int A_rows, unsigned int lda,const float *B, unsigned int B_cols, unsigned int ldb, float *d_C, unsigned int ldc, const float alpha = 1.0f, const float beta = 1.0f);

