#pragma once
#include "comp_mats.h"
#include <string>
#include <map>
#include <set>
#include <vector>


typedef std::map<intT, std::set<intT> > GraphMap;

typedef std::vector<intT> svi;

typedef std::vector<DataT> svd;


//permutes an array of n elements (original) according to a permutation (perm);
//es 
// array: A B C D E
// perm: 
template <class myType>
int permute(myType* arr, intT* perm, intT n) {
	intT i = 0;
	myType temp;
	intT cur_row, next_row;
	//Do the permutations cycle by cycle, until all elements have been permuted.
	while (i < n) {//LI: all elements up to i appeared in a cycle
		temp = arr[i];
		intT start_i = i;
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
		while ((i < n) and (perm[i] < 0)) i++;
	}

	//reconstruct perm (remove flags)
	for (i = 0; i < n; i++) {
		perm[i] = -perm[i] - 1;
	}
	return 0;

}

template <class myType>
myType IDX(myType row, myType col, myType lead_dim, int fmt)
{
	//finds storing index of a matrix elements given
	//row: row index
	//col: columnn index
	//lead_dimension: the leading storing dimension
	//fmt:      0: row major
	//          1: column major

	if (fmt == 0)
	{
		return row * lead_dim + col;
	}
	else
	{
		return col * lead_dim + row;
	}
}

template <class myType>
myType leading_dim(myType rows, myType cols, int fmt)
{
	return (fmt == 0) ? cols : rows;
}

int is_empty(DataT* mat, intT rows, intT cols, intT lead_dim, int fmt);

int mat_cpy(DataT* in_mat, intT in_rows, intT in_cols, intT in_lead_dim, int in_fmt, DataT* out_mat, intT out_lead_dim, int out_fmt);

int random_mat(DataT* mat, intT rows, intT cols, float sparsity);

int equal(intT rows, intT cols, DataT* A, intT lead_A, int fmt_A, DataT* B, intT lead_B, int fmt_B, DataT eps);

int random_sparse_blocks_mat(DataT* mat, intT rows, intT cols, int fmt, intT block_size, float block_sparsity, float block_entries_sparsity);

int random_sparse_blocks_mat(VBS& vbmat, intT rows, intT cols, int blocks_fmt, int entries_fmt, intT row_block_size, intT col_block_size, float block_density, float entries_density);

int matprint(DataT* mat, intT rows, intT cols, intT lead_dim, int fmt, bool struct_only);

int matprint(DataT* mat, intT rows, intT* row_part, intT row_blocks, intT cols, intT* col_part, intT col_blocks, intT lead_dim, int fmt, bool struct_only);

int arr_print(intT* arr, intT len);

int sort_permutation(intT* perm, intT* arr, intT n);

int partition(intT* arr, intT start, intT end, int step);

int randperm(intT* arr, intT len);

intT* rand_partition(intT* part, intT len, intT blocks);

int cleanVBS(VBS& vbmat);

int init_VBS(VBS& vbmat, intT block_rows, intT* row_part, intT block_cols, intT* col_part, int blocks_fmt, int entries_fmt);

int convert_to_VBS(DataT* mat, intT mat_rows, intT mat_cols, int mat_fmt, VBS& vbmat, intT block_rows, intT* row_part, intT block_cols, intT* col_part, int vbmat_blocks_fmt, int vbmat_entries_fmt, int no_zero_mode);

int convert_to_mat(const VBS& vbmat, DataT* out_mat, int out_mat_fmt);

int convert_to_VBS(const CSR& cmat, VBS& vbmat, intT block_rows, intT* rowpart, intT block_cols, intT* colpart, int vbmat_block_fmt, int vbmat_entries_fmt);

int convert_to_VBS(const CSR& cmat, VBS& vbmat, intT row_block_size, intT col_block_size, int vbmat_block_fmt, int vbmat_entries_fmt);

int matprint(const VBS& vbmat);

int cleanCSR(CSR& cmat);

int convert_to_mat(const CSR& cmat, DataT* out_mat, int out_mat_fmt);

int convert_to_CSR(const DataT* in_mat, intT mat_rows, intT mat_cols, int mat_fmt, CSR& cmat, int cmat_fmt);

int convert_to_CSR(const VBS& vbmat, CSR& cmat, int csr_fmt);

int matprint(const CSR& cmat);

int copy(const CSR& in_cmat, CSR& out_cmat);

int transpose(const CSR& in_cmat, CSR& out_cmat, int new_fmt);

int permute_CSR(CSR& cmat, intT* perm, int dim);

intT count_nnz(CSR& cmat);

intT count_nnz_blocks(VBS& vbmat);

intT scalar_product(intT* pat_0, intT len_0, intT* pat_1);

intT norm2(intT* arr, intT len);

void read_snap_format(GraphMap& gmap, std::string filename);

void read_snap_format(GraphMap& gmap, std::string filename, std::string delimiter);

int isProper(const GraphMap& gmap, bool mirroring);

void MakeUndirected(GraphMap& gmap);

void MakeProper(GraphMap& gmap);

void write_snap_format(GraphMap& gmap, std::string filename);

void convert_to_CSR(const GraphMap& gmap, CSR& cmat, int cmat_fmt);

void read_edgelist(std::string filename, CSR& cmat, int cmat_fmt, std::string delimiter);

