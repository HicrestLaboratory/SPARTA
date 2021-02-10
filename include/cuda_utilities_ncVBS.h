#pragma once

void cublas_ncVBS_multiply(ncVBS& vbmatA, const DataT* B, int B_cols, int B_lead_dim, DataT* C, int C_lead_dim, float& dt, int n_streams_mult, int n_streams_cpy);
