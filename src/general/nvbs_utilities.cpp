#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream> //ifstream

#include <cmath> // std::abs
#include <string>
#include <map>
#include <set>
#include <string>

#include <random>
#include <vector>
#include <algorithm>    // std::random_shuffle

#include "sparse_utilities.h"
#include "nvbs_utilities.h"

#include "comp_mats.h"


//REMOVE!!!!!!!
typedef int intT;
typedef float DataT;
//!!!!!!!!!!!!!

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
    vbmat.nzcount = new intT[block_cols]{ 0 };
    vbmat.nzindex = new intT * [block_cols];
    vbmat.mab = new DataT * [block_cols];
}



//INSPECTION UTILITIES

int matprint(const ncVBS& vbmat)
{
    for (intT jb = 0; jb < vbmat.block_cols; jb++)
    {
        intT rows_number = vbmat.nzcount[jb];
        intT* rows_indices = vbmat.nzindex[jb];
        intT column_start = vbmat.col_part[jb];
        intT column_end = vbmat.col_part[jb + 1];
        intT column_block_size = column_end - column_start;

        std::cout << "BLOCK " << jb << "; rows: " << rows_number << "; start: " << column_start << "; end: " << column_end << std::endl;

        matprint(vbmat.mab[jb], rows_number, column_block_size, column_block_size, 0);
        
        std::cout << "INDICES" << std::endl;
        arr_print(vbmat.nzindex[jb], rows_number);

    }



}

int print_block_info(const ncVBS& vbmat, int block)
{
    //NOT IMPLEMENTED YET
    return 0;
}

bool equal(ncVBS& vbmat, DataT* mat, intT mat_rows, intT mat_cols, intT mat_leading_dim, int mat_fmt)
{
    DataT* mat2 = new DataT[vbmat.rows * vbmat.cols()];
    convert_to_mat(vbmat, mat2, 0);
    return equal(vbmat.rows, vbmat.cols(), mat2, vbmat.cols(), 0, mat, mat_leading_dim, mat_fmt, 0.00001f);
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
    intT nonzero_rows_per_block = (int)(row_density * mat_rows);
    intT elements_per_row = (int)(mat_density * block_size);
    intT elements_per_block = nonzero_rows_per_block * block_size;
    if (mat_cols % block_size != 0) std::cout << "WARNING: mat_cols not a multiple of block_size" << std::endl;


    vbmat.rows = mat_rows;
    intT block_cols = mat_cols / block_size;

    
    vbmat.block_cols = block_cols;
    vbmat.col_part = new intT[block_cols + 1];
    partition(vbmat.col_part, 0, mat_cols, block_size);

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
            if (rows[i] == 1)
            {
                vbmat.nzindex[jb][nz_row] = i;
                nz_row++;
            }
        }

        arr_print(vbmat.nzindex[jb], nonzero_rows_per_block);

        //fill mab[jb] with random rows
        for (intT i = 0; i < nonzero_rows_per_block; i++)
        {

            //add one random row to mab[jb]
            svi elems = svi(block_size, 0);
            std::fill(elems.begin(), elems.begin() + elements_per_row, 1); //only elems_per_row are nonzero;
            std:random_shuffle(elems.begin(), elems.end());
            std::copy(elems.begin(), elems.end(), vbmat.mab[jb] + (block_size * i)); 
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

int convert_to_mat(ncVBS& vbmat, DataT* out_mat, int C_fmt)
{
    //C_mat must be of the appropriate dimension; 

    intT out_mat_rows = vbmat.rows;
    intT out_mat_cols = vbmat.cols();
    intT mat_leading_dim = C_fmt == 0 ? out_mat_cols : out_mat_rows;

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
                intT mat_idx = IDX(i, j, mat_leading_dim, C_fmt);

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
    intT mat_rows = cmat.rows;
    intT mat_cols = cmat.cols;

    if (cmat.fmt != 0)
    {
        std::cout << "ERROR! CSR must be in row-major fmt (fmt = 0)" << std::endl;
    }

    initialize_ncVBS(vbmat, mat_rows, block_cols, col_part);

    std::vector<intT> nzindices[block_cols];
    std::vector<DataT> mab[block_cols];



    //count nonzero rows per block;
    for (intT i = 0; i < mat_rows; i++)
    {
        intT block = 0;
        intT n_elems = cmat.nzcount[i];
        intT nz = 0;
        while(nz < n_elems)
        {
            DataT elem = cmat.ma[i][nz];
            while (cmat.ja[i][nz] >= col_part[block + 1]) block++; //find in which block the nz is;
            vbmat.nzcount[block]++; //flag this row as nonzero;

            intT width = vbmat.block_width(block);
            nzindices[block].push_back(i); //store the row index;
            std::cout << "i: " << i << " width: " << width << "nz: " << cmat.ja[i][nz] << " block: " << block << std::endl;
            std::cout << "start " << col_part[block] << " end " << col_part[block + 1] << std::endl;
            std::vector<DataT> temp_vec(width);

            while (nz < n_elems && cmat.ja[i][nz] < col_part[block + 1])
            {
                intT column = cmat.ja[i][nz] - col_part[block];
                DataT elem = cmat.ma[i][nz];
                temp_vec[column] = elem;

                nz++;
            }
            mab[block].insert(mab[block].end(), temp_vec.begin(), temp_vec.end());
        }
    }


    //initialize vbmat.nzindex;
    for (intT jb = 0; jb < block_cols; jb++)
    {
        intT width = vbmat.block_width(jb);
        std::cout << "jb: " << jb << "part start: " << col_part[jb] << " part end: " << col_part[jb+1] << " width: "<< width << "mab len: " << mab[jb].size() << std::endl;

        vbmat.nzindex[jb] = new intT[vbmat.nzcount[jb]]{ 0 };
        vbmat.mab[jb] = new DataT[vbmat.nzcount[jb] * width]{ 0 };
        std::copy(nzindices[jb].begin(), nzindices[jb].end(), vbmat.nzindex[jb]);
        std::copy(mab[jb].begin(), mab[jb].end(), vbmat.mab[jb]);
    }

}

int convert_to_CSR(const ncVBS& vbmat, CSR& cmat, int csr_fmt)
{
    //NOT IMPLEMENTED YET
    return 0;
}




//MULTIPLICATION UTILITY

int multiply(const ncVBS& vbmat, DataT* B_mat, intT B_mat_cols, int B_mat_fmt, intT B_leading_dim, DataT* C_mat, intT C_leading_dim, int C_fmt)
{

    for (intT jb = 0; jb < vbmat.block_cols; jb++)
    {
        intT rows_number = vbmat.nzcount[jb];
        intT* rows_indices = vbmat.nzindex[jb];
        intT column_start = vbmat.col_part[jb];
        intT column_end = vbmat.col_part[jb + 1];
        intT column_block_size = column_end - column_start;

        for (int nz_i = 0; nz_i < rows_number; nz_i++)
        {
            int i = rows_indices[nz_i];
            for (int col_B = 0; col_B < B_mat_cols; col_B++)
            {
                DataT elem = 0;
                intT C_IDX = IDX(i, col_B, C_leading_dim, C_fmt);
                for (int j = column_start; j < column_end; j++)
                {
                    intT B_IDX = IDX(j, col_B, B_leading_dim, B_mat_fmt);
                    intT vbmat_IDX = IDX(nz_i, j - column_start, column_block_size, 0);
                    elem += B_mat[B_IDX] * vbmat.mab[jb][vbmat_IDX];
                }
                C_mat[C_IDX] += elem;
            }
        }

    }
    return 0;
}

int multiply(DataT* A_mat, intT A_mat_rows, intT A_mat_cols, int A_mat_fmt, intT A_mat_leading_dim, DataT* B_mat, intT B_mat_cols, int B_mat_fmt, intT B_mat_leading_dim, DataT* C_mat, intT C_mat_leading_dim, int C_mat_fmt)
{
    for (intT i = 0; i < A_mat_rows; i++)
    {
        for (intT col_B = 0; col_B < B_mat_cols; col_B++)
        {
            DataT elem = 0;
            intT C_IDX = IDX(i, col_B, C_mat_leading_dim, C_mat_fmt);
            for (intT j = 0; j < A_mat_cols; j++)
            {
                intT A_IDX = IDX(i, j, A_mat_leading_dim, A_mat_fmt);
                intT B_IDX = IDX(j, col_B, B_mat_leading_dim, B_mat_fmt);
                elem += A_mat[A_IDX] * B_mat[B_IDX];
            }
            C_mat[C_IDX] = elem;
        }
    }
}