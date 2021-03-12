#pragma once
#include <string>
#include "comp_mats.h"

void cublas_ncVBS_multiply(ncVBS& vbmatA, const DataT* B, int B_cols, int B_lead_dim, DataT* C, int C_lead_dim, float* times, int n_streams_mult = 32, int n_streams_cpy = 32);
