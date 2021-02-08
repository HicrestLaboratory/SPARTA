#pragma once
#include "comp_mats.h"
#include "sparse_utilities.h"

int cleanVBS(ncVBS& vbmat);

int initialize_ncVBS(ncVBS& vbmat, intT mat_rows, intT block_cols, intT* col_part);

int random_ncVBS(ncVBS& vbmat, intT mat_rows, intT mat_cols, intT block_size, float mat_density, float row_density, int mode = 0, int block_variation = 0, float density_variation = 0);

int convert_to_mat(ncVBS& vbmat, DataT* out_mat, int out_mat_fmt);

int convert_to_ncVBS(DataT* mat, intT mat_rows, intT mat_cols, int mat_fmt, int mat_lead_dim, ncVBS& vbmat, intT block_cols, intT* col_part);

int matprint(const ncVBS& vbmat);

bool equal(ncVBS& vbmat, const DataT* mat, intT mat_rows, intT mat_cols, intT mat_leading_dim, int mat_fmt);

int multiply(const ncVBS& vbmat, DataT* in_mat, intT in_mat_cols, int in_mat_fmt, intT in_mat_leading_dim, DataT* out_mat, intT out_mat_leading_dim, int out_mat_fmt);

int multiply(DataT* A_mat, intT A_mat_rows, intT A_mat_cols, int A_mat_fmt, intT A_mat_leading_dim, DataT* B_mat, intT B_mat_cols, int B_mat_fmt, intT B_mat_leading_dim, DataT* C_mat, intT C_mat_leading_dim, int C_mat_fmt);
