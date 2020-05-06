#pragma once

typedef float DataT;

int IDX(int row, int col, int lead_dim, int fmt);

int is_empty(DataT* mat, int rows, int cols, int lead_dim, int fmt);

int mat_cpy(DataT* in_mat, int in_rows, int in_cols, int in_lead_dim, int in_fmt, DataT* out_mat, int out_lead_dim, int out_fmt);

int random_mat(DataT* mat, int rows, int cols, float sparsity);

int random_sparse_blocks_mat(DataT* mat, int rows, int cols, int fmt, int block_size, float block_sparsity, float block_entries_sparsity);

int matprint(DataT* mat, int rows, int cols, int lead_dim, int fmt);

int cleanVBS(VBSfx& vbmat);

int convert_to_VBSfx(DataT* mat, int mat_rows, int mat_cols, int mat_fmt, VBSfx& vbmat, int block_size, int vbmat_block_fmt, int vbmat_entries_fmt);

int convert_to_mat(const VBSfx& vbmat, DataT* out_mat, int out_mat_fmt);

int convert_to_VBSfx(const CSR& cmat, VBSfx& vbmat, int block_size, int vbmat_block_fmt, int vbmat_entries_fmt);

int cleanCSR(CSR& cmat);

int convert_to_mat(const CSR& cmat, DataT* out_mat, int out_mat_fmt);

int convert_to_CSR(const DataT* in_mat, int mat_rows, int mat_cols, int mat_fmt, CSR& cmat, int csr_fmt);

int convert_to_CSR(const VBSfx& vbmat, CSR& cmat, int csr_fmt);

int transpose(const CSR& in_cmat, CSR& out_cmat, int fmt_change);

int hash_permute(CSR& cmat, int block_size, int* perm, int* group, int mode);

int hash(int* arr, int a_len, int block_size, int mode);

int check_same_pattern(int* arr0, int len0, int* arr1, int len1, int block_size, int mode);

int permute(CSR& cmat, int* perm, int dim);
