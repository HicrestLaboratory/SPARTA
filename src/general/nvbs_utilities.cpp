#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream> //ifstream

#include <cmath> // std::abs
#include <string>
#include <map>
#include <set>

#include <random>
#include <vector>
#include <algorithm>    // std::random_shuffle

#include "sparse_utilities.h"
#include "nvbs_utilities.h"

#include "comp_mats.h"

int cleanVBS(ncVBS& vbmat)
{

    if (vbmat.nzcount) delete[] vbmat.nzcount;

    if (vbmat.mab)
    {
        for (int i = 0; i < block_cols; i++)
        {
            if (vbmat.mab[i]) delete[] vbmat.mab[i];
        }
        delete[] vbmat.mab;
    }

    if (vbmat.nzindex)
    {
        for (int i = 0; i < block_cols; i++)
        {
            if (vbmat.nzindex[i]) delete[] vbmat.nzindex[i];
        }
        delete[] vbmat.nzindex;
    }

    if (vbmat.col_part) delete[] vbmat.col_part;

    return 0;
}



//INSPECTION UTILITIES

int matprint(const ncVBS& vbmat)
{
    //NOT IMPLEMENTED YET
    return 0;
}

int print_block_info(const ncVBS& vbmat, int block)
{
    //NOT IMPLEMENTED YET
    return 0;
}

bool equal(const ncVBS& vbmat_1, const ncVBS& vbmat_2)
{
    //NOT IMPLEMENTED YET
    return false;
}

bool equal(const ncVBS& vbmat, const DataT* mat, intT mat_rows, intT mat_cols, intT mat_leading_dim, int mat_fmt)
{
    //NOT IMPLEMENTED YET
    return false;
}

bool check_consistent(const ncVBS& vbmat)
{
    //NOT IMPLEMENTED YET
    return false;
}

int count_nnz(const ncVBS& vbmat)
{
    //NOT IMPLEMENTED YET
    return -1;
}

int count_nnz_rows(const ncVBS& vbmat)
{
    //NOT IMPLEMENTED YET
    return -1;
}




//GENERATION UTILITIES

int random_ncVBS(ncVBS& vbmat, intT mat_rows, intT mat_cols, intT block_size, float mat_density, float row_density, int mode = 0, int block_variation = 0, float density_variation = 0)
{
    //NOT IMPLEMENTED YET
    return 0;
}

int random_ncVBS_partitioned(ncVBS& vbmat, intT mat_rows, intT mat_cols, intT block_cols, intT* col_part, float mat_density, float row_density, int mode = 0, float density_variation = 0)
{
    //NOT IMPLEMENTED YET
    return 0;
}



//CONVERSION UTILITIES

int convert_to_mat(const ncVBS& vbmat, DataT* out_mat, int out_mat_fmt)
{
    //out_mat must be of the appropriate dimension; 

    intT out_mat_rows = vbmat.rows;
    intT out_mat_cols = vbmat.cols();
    mat_leading_dim = out_mat_fmt == 0 ? out_mat_cols : out_mat_rows;

    for (int jb = 0; jb < vbmat.block_cols; jb++)
    {
        intT rows_number = vbmat.nzcount[i];
        intT* rows_indices = vbmat.nzindex[i];
        intT column_start = vbmat.col_part[i];
        intT column_end = vbmat.col_part[i + 1];
        intT column_block_size = column_end - column_start;

        for (int row = 0; row < rows_number; row++)
        {
            intT i = rows_indices[row];
            for (int j = column_start; j < column_end; j++)
            {
                intT vbmat_idx = IDX(row, j - column_start, column_block_size, 0);
                intT mat_idx = IDX(i, j, mat_leading_dim, out_mat_fmt);

                out_mat[mat_idx] = vbmat.mab[vbmat_idx];
            }
        }

    }
    return 0;
}

int convert_to_ncVBS(DataT* mat, intT mat_rows, intT mat_cols, int mat_fmt, ncVBS& vbmat, intT block_cols, intT* col_part)
{
    //NOT IMPLEMENTED YET
    return 0;
}

int convert_to_ncVBS(const CSR& cmat, ncVBS& vbmat, intT block_cols, intT* col_part)
{
    //NOT IMPLEMENTED YET
    return 0;
}

int convert_to_CSR(const ncVBS& vbmat, CSR& cmat, int csr_fmt)
{
    //NOT IMPLEMENTED YET
    return 0;
}



//MULTIPLICATION UTILITY

int multiply(const ncVBS& vbmat, DataT* in_mat, intT in_mat_rows, intT in_mat_cols, int in_mat_fmt, intT in_mat_leading_dim, DataT* out_mat, intT out_mat_leading_dim, int out_mat_fmt)
{
    //NOT IMPLEMENTED YET
    return 0;
}