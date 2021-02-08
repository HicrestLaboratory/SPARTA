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
        for (int i = 0; i < vbmat.block_cols; i++)
        {
            if (vbmat.mab[i]) delete[] vbmat.mab[i];
        }
        delete[] vbmat.mab;
    }

    if (vbmat.nzindex)
    {
        for (int i = 0; i < vbmat.block_cols; i++)
        {
            if (vbmat.nzindex[i]) delete[] vbmat.nzindex[i];
        }
        delete[] vbmat.nzindex;
    }

    if (vbmat.col_part) delete[] vbmat.col_part;

    return 0;
}
int initialize_ncVBS(ncVBS& vbmat, intT mat_rows, intT block_cols, intT* col_part)
{
    vbmat.block_cols = block_cols;
    vbmat.col_part = new intT[block_cols + 1];
    std::copy(col_part, col_part + block_cols + 1, vbmat.col_part);
    vbmat.rows = mat_rows;
    vbmat.nzcount = new intT[block_cols];
    vbmat.nzindex = new intT * [block_cols];
    vbmat.mab = new DataT * [block_cols];
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
    intT nonzero_rows_per_block = (int)row_density * mat_rows;
    intT elements_per_row = (int)mat_density * block_size;
    intT elements_per_block = elements_per_row * nonzero_rows_per_block;
    if (mat_cols % block_size != 0) std::cout << "WARNING: mat_cols not a multiple of block_size" << std::endl;

    vbmat.rows = mat_rows;
    intT block_cols = mat_cols / block_size;
    vbmat.block_cols = block_cols;
    vbmat.col_part = new intT[block_cols + 1];
    partition(vbmat.col_part, 0, block_cols, block_size);

    vbmat.nzcount = new intT[block_cols];
    vbmat.nzindex = new intT * [block_cols];
    vbmat.mab = new DataT * [block_cols];

    //loop through column blocks
    for (int jb = 0; jb < block_cols; jb++)
    {
        vbmat.nzcount[jb] = nonzero_rows_per_block;
        vbmat.nzindex[jb] = new intT[vbmat.nzcount[jb]];
        vbmat.mab[jb] = new DataT[elements_per_block];

        //determine nonzero rows for a block
        svi rows = svi(mat_rows, 0);              //will store 0 unless a row is nonzero;
        std::fill(rows.begin(), rows.begin() + nonzero_rows_per_block, 1);    //make nonzero_rows_per_block rows nonzero;
        std::random_shuffle(rows.begin(), rows.end());          //put the nonzero rows in random positions

        //fill nzindex[jb] with random indices
        intT nz_row = 0;
        for (intT i = 0; i < mat_rows; i++)
        {
            if (rows[i] == 1) vbmat.nzindex[jb][nz_row] = i;
            nz_row++;
        }

        //fill mab[jb] with random rows
        for (intT i = 0; i < nonzero_rows_per_block; i++)
        {
            //add one random row to mab[jb]
            svi elems = svi(block_size, 0);
            std::fill(elems.begin(), elems.begin() + elements_per_row, 1); //only elems_per_row are nonzero;
            std:random_shuffle(rows.begin(), rows.end());
            std::copy(rows.begin(), rows.end(), vbmat.mab[jb] + block_size * i); 
        }

    }

    return 0;
}

int random_ncVBS_partitioned(ncVBS& vbmat, intT mat_rows, intT mat_cols, intT block_cols, intT* col_part, float mat_density, float row_density, int mode = 0, float density_variation = 0)
{
    //NOT IMPLEMENTED YET
    return 0;
}



//CONVERSION UTILITIES

int convert_to_mat(ncVBS& vbmat, DataT* out_mat, int out_mat_fmt)
{
    //out_mat must be of the appropriate dimension; 

    intT out_mat_rows = vbmat.rows;
    intT out_mat_cols = vbmat.cols();
    intT mat_leading_dim = out_mat_fmt == 0 ? out_mat_cols : out_mat_rows;

    std::cout << vbmat.block_cols << std::endl;
    for (int jb = 0; jb < vbmat.block_cols; jb++)
    {
        intT rows_number = vbmat.nzcount[jb];
        intT* rows_indices = vbmat.nzindex[jb];
        intT column_start = vbmat.col_part[jb];
        intT column_end = vbmat.col_part[jb + 1];
        intT column_block_size = column_end - column_start;

        for (int row = 0; row < rows_number; row++)
        {
            intT i = rows_indices[row];
            for (int j = column_start; j < column_end; j++)
            {
                intT vbmat_idx = IDX(row, j - column_start, column_block_size, 0);
                intT mat_idx = IDX(i, j, mat_leading_dim, out_mat_fmt);

                out_mat[mat_idx] = vbmat.mab[jb][vbmat_idx];
            }
        }

    }
    return 0;
}

int convert_to_ncVBS(DataT* mat, intT mat_rows, intT mat_cols, int mat_fmt, int mat_lead_dim, ncVBS& vbmat, intT block_cols, intT* col_part)
{
    initialize_ncVBS(vbmat, mat_rows, block_cols, col_part);

    for (int jb = 0; jb < block_cols; jb++)
    {
        intT block_start = col_part[jb];
        intT block_end = col_part[jb + 1];
        intT block_size = block_end - block_start;

        svi nz_rows; 
        //fill nz_rows with indices of nonzero rows for this block;
        for (intT i = 0; i < mat_rows; i++)
        {

            //check if row is empty
            bool empty = true;
            for (int j = block_start; j < block_end; j++)
            {
                intT mat_idx = IDX(i, j, mat_lead_dim, mat_fmt);
                if (mat[mat_idx] != 0)
                {
                    empty = false;
                    break;
                }
            }

            if (!empty) nz_rows.push_back(i);
        }


        //initialize and fill nzindex[jb];
        intT row_size = nz_rows.size();
        vbmat.nzcount[jb] = row_size;
        if (row_size > 0)
        {
            vbmat.nzindex[jb] = new intT[row_size];
            std::copy(nz_rows.begin(), nz_rows.end(), vbmat.nzindex[jb]);

        }


        //initialize and fill mab[jb]
        if (row_size > 0)
        {
            intT number_of_elements = row_size * block_size;
            vbmat.mab[jb] = new DataT[number_of_elements];
    
            for (int i = 0; i < row_size; i++)
            {
                int i_mat = nz_rows[i];
                for (int j = 0; j < block_size; j++)
                {
                    int j_mat = j + block_start;
                    int mat_idx = IDX(i_mat, j_mat, mat_lead_dim, mat_fmt);
                    int vbmat_idx = IDX(i, j, block_size, 0);
                    vbmat.mab[jb][vbmat_idx] = mat[mat_idx];
                }
            }
        }

    }

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