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

int cleanVBS(VBS& vbmat);

int convert_to_VBS(DataT* mat, int mat_rows, int mat_cols, int mat_fmt, VBS& vbmat, int block_rows, int* rowpart, int block_cols, int* colpart, int vbmat_blocks_fmt, int vbmat_entries_fmt);

int convert_to_mat(const VBS& vbmat, DataT* out_mat, int out_mat_fmt);

int convert_to_VBS(const CSR& cmat, VBS& vbmat, int block_rows, int* rowpart, int block_cols, int* colpart, int vbmat_block_fmt, int vbmat_entries_fmt);

int matprint(const VBS& vbmat);

int cleanCSR(CSR& cmat);

int convert_to_mat(const CSR& cmat, DataT* out_mat, int out_mat_fmt);

int convert_to_CSR(const DataT* in_mat, int mat_rows, int mat_cols, int mat_fmt, CSR& cmat, int csr_fmt);

int convert_to_CSR(const VBSfx& vbmat, CSR& cmat, int csr_fmt);

int convert_to_CSR(const VBS& vbmat, CSR& cmat, int csr_fmt);

int transpose(const CSR& in_cmat, CSR& out_cmat, int fmt_change);

int hash(int* arr, int a_len, int block_size, int mode);

int hash(int* arr, int a_len, int* block_partition, int mode);

int check_same_pattern(int* arr0, int len0, int* arr1, int len1, int block_size, int mode);

int check_same_pattern(int* arr0, int len0, int* arr1, int len1, int* block_partition, int mode);

int pattern_scalar_product(int* arr0, int len0, int* pattern, int* block_partition, int mode);

int permute_CSR(CSR& cmat, int* perm, int dim);

int hash_permute(CSR& cmat, int* comp_dim_partition, int* perm, int* group, int mode);

//permutes an array of n elements (original) according to a permutation (perm);
template <class myType>
int permute(myType* arr, int* perm, int n) {
	int i = 0;
	myType temp;
	int cur_row, next_row;

	//Do the permutations cycle by cycle, until all elements have been permuted.
	while (i < n) {//LI: all elements up to i appeared in a cycle
		temp = arr[i];
		int start_i = i;
		next_row = i;

		//one cycle of permutations
		while (perm[next_row] >= 0) {//check if we are back to the start of the cycle
			cur_row = next_row;
			next_row = perm[cur_row];
			arr[cur_row] = arr[next_row];
			perm[cur_row] = (perm[cur_row] + 1) * (-1);//flag perm element as executed 
		}
		//end cycle and procede to next element
		if (next_row != start_i) return 1; //PERM IS NOT A PROPER PERMUTATION
		arr[cur_row] = temp;
		while (perm[i] < 0) i++;
	}

	//reconstruct perm (remove flags)
	for (i = 0; i < n; i++) {
		perm[i] = -perm[i] - 1;
	}
	return 0;

}

