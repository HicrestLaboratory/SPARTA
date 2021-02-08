#pragma once
#include "comp_mats.h"


int cleanVBS(ncVBS& vbmat);

int initialize_ncVBS(ncVBS& vbmat, intT mat_rows, intT block_cols, intT* col_part);

int random_ncVBS(ncVBS& vbmat, intT mat_rows, intT mat_cols, intT block_size, float mat_density, float row_density, int mode, int block_variation, float density_variation);

int convert_to_mat(const ncVBS& vbmat, DataT* out_mat, int out_mat_fmt);

int convert_to_ncVBS(DataT* mat, intT mat_rows, intT mat_cols, int mat_fmt, int mat_lead_dim, ncVBS& vbmat, intT block_cols, intT* col_part);
