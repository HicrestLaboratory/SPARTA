#pragma once
#include "comp_mats.h"


//permutes an array of n elements (original) according to a permutation (perm);
//es 
// array: A B C D E
// perm: 

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

int IDX(int row, int col, int lead_dim, int fmt);

int is_empty(DataT* mat, int rows, int cols, int lead_dim, int fmt);

int mat_cpy(DataT* in_mat, int in_rows, int in_cols, int in_lead_dim, int in_fmt, DataT* out_mat, int out_lead_dim, int out_fmt);

int random_mat(DataT* mat, int rows, int cols, float sparsity);

int leading_dim(int rows, int cols, int fmt);

int equal(int rows, int cols, DataT* A, int lead_A, int fmt_A, DataT* B, int lead_B, int fmt_B, DataT eps);

int random_sparse_blocks_mat(DataT* mat, int rows, int cols, int fmt, int block_size, float block_sparsity, float block_entries_sparsity);

int matprint(DataT* mat, int rows, int cols, int lead_dim, int fmt);

int matprint(DataT* mat, int rows, int* row_part, int row_blocks, int cols, int* col_part, int col_blocks, int lead_dim, int fmt);

int arr_print(int* arr, int len);

int sort_permutation(int* perm, int* arr, int n);

int* linspan(int start, int end, int step);

int* randperm(int len);

int* rand_partition(int* part, int len, int blocks);

int count_groups(int* grp, int grp_len);

int grp_to_partition(int* grp, int grp_len, int* partition);

int cleanVBS(VBS& vbmat);

int convert_to_VBS(DataT* mat, int mat_rows, int mat_cols, int mat_fmt, VBS& vbmat, int block_rows, int* row_part, int block_cols, int* col_part, int vbmat_blocks_fmt, int vbmat_entries_fmt, int no_zero_mode = 0);

int convert_to_mat(const VBS& vbmat, DataT* out_mat, int out_mat_fmt);

int convert_to_VBS(const CSR& cmat, VBS& vbmat, int block_rows, int* rowpart, int block_cols, int* colpart, int vbmat_block_fmt, int vbmat_entries_fmt, int no_zero_mode = 0);

int matprint(const VBS& vbmat);

int cleanCSR(CSR& cmat);

int convert_to_mat(const CSR& cmat, DataT* out_mat, int out_mat_fmt);

int convert_to_CSR(const DataT* in_mat, int mat_rows, int mat_cols, int mat_fmt, CSR& cmat, int cmat_fmt);

int convert_to_CSR(const VBS& vbmat, CSR& cmat, int csr_fmt);

int matprint(const CSR& cmat);

int copy(const CSR& in_cmat, CSR& out_cmat);

int transpose(const CSR& in_cmat, CSR& out_cmat, int new_fmt);

int permute_CSR(CSR& cmat, int* perm, int dim);

int count_nnz(CSR& cmat);

int hash_permute(CSR& cmat, int* comp_dim_partition, int* perm, int* group, int mode);

int hash(int* arr, int a_len, int block_size, int mode);

int hash(int* arr, int a_len, int* block_partition, int mode);

int check_same_pattern(int* arr0, int len0, int* arr1, int len1, int block_size, int mode);

int check_same_pattern(int* arr0, int len0, int* arr1, int len1, int* block_partition, int mode);

int get_pattern(int* arr0, int len0, int* block_partition, int* pattern, int mode);

int scalar_product(int* pat_0, int len_0, int* pat_1);

int norm2(int* arr, int len);

int angle_method(CSR& cmat, float eps, int* comp_dim_partition, int nB, int* in_perm, int* in_group, int* out_group, int mode);
