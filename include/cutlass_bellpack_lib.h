#pragma once

#define CUTLASS
#ifdef CUTLASS

#include "definitions.h"

template<typename iT, typename T>
int compute_cutlass_bellpack (int rows, int cols, int ell_blocksize, int ellValue_cols, iT* ellColInd, T *ellValues, int B_rows, int B_cols, T *B_vals, int C_rows, int C_cols, T *C_vals);

void bellpack_cutlass_multiplyAB(VBR* A, DataT* B, int B_cols, DataT_C* C, int C_cols, float& dt, int verbose);

#endif
