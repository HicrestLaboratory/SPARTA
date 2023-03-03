#pragma once

#include "matrices.h"

#define PICO_DEBUG
// #define B_ROW_GROWING
// #define B_COL_GROWING
#define BDG_CKP { printf("BDG_CKP: file %s line %d\n", __FILE__, __LINE__); }

void cublas_blockmat_multiply(const VBR& vbmatA, DataT* B, int B_rows, int B_lead_dim, DataT_C* C, int C_lead_dim, float &dt, int n_streams = 16);

void cublas_blockmat_multiplyAB(const VBR& vbmatA, DataT* B, int B_cols, DataT_C* C, float& dt, int n_streams = 16);

//int cublas_gemm_custom(const DataT* A, unsigned int A_rows, unsigned int A_cols, unsigned int lda, const DataT* B, unsigned int B_cols, unsigned int ldb, DataT_C* C, unsigned int ldc, const DataT_C alpha, const DataT_C beta, float& dt);

int cusparse_gemm_custom(int rows, int cols, int nnz, int* csrRowPtr, int* csrColInd, DataT* csrVal, DataT* B, int B_cols, int B_lead_dim, DataT_C* C, int C_lead_dim, const DataT_C alpha, const DataT_C beta, float& dt);

void pico_print_SpMMM(const char* Aname, int An, int Am, int Az, int* Arows, int* Acols, DataT* Avals, const char* Bname, int Bn, int Bm, DataT* B, const char* Cname, long int Cn, long int Cm, DataT_C* C);

void pico_print_SpMMM(const char* Aname, VBR* A, const char* Bname, int Bn, int Bm, DataT* B, const char* Cname, long int Cn, long int Cm, DataT_C* C);

void pico_print_SpMMM(const char* Aname, int rows, int cols, int ell_blocksize, int ellValue_cols, int ellColumnsInd_rows, int ellColumnsInd_cols, int num_blocks, intT* ellColumnsInd, DataT_C* ellValues, const char* Bname, int Bn, int Bm, DataT* B, const char* Cname, long int Cn, long int Cm, DataT_C* C);

int cusparse_gemm_custom2(int rows, int cols, int nnz, int* csrRowPtr, int* csrColInd, DataT* csrVal, DataT* B, int B_cols, int B_lead_dim, DataT_C* C, int C_lead_dim, const DataT_C alpha, const DataT_C beta, float& dt);

int prepare_cusparse_CSR(CSR& cmat, int **csrRowPtr, int **csrColInd, DataT **csrVal);

int cusparse_gemm_custom_ellpack(int rows, int cols, int A_ell_blocksize, int A_ellValues_cols, int A_ellColInd_cols, int A_ellColInd_rows, int A_num_blocks, intT* A_ellColInd, DataT_C* A_ellValues, DataT* B, int B_cols, int B_lead_dim, DataT_C* C, int C_lead_dim, const DataT_C alpha, const DataT_C beta, float& dt);

int prepare_cusparse_BLOCKEDELLPACK(VBR *A, int *ell_blocksize, int* ellValue_cols, int *ellColInd_rows, int *ellColInd_cols, int *num_blocks, intT** ellColInd, DataT_C** ellValues);

