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
#include "comp_mats.h"
#include "omp.h"


// Matrix utilities

std::string print_mat_val(DataT val, bool struct_only, std::string zero_style, std::string nz_style)
{
    std::string out;
    if (struct_only)
    {
        if (val == 0) out = zero_style;
        else out = nz_style;
    }
    else
    {
        out = std::to_string(val);
    }
    return out;
}

int is_empty(DataT* mat, intT rows, intT cols, intT lead_dim, int fmt) {
    //check if a matrix is empty

    for (intT i = 0; i < rows; i++)
    {
        for (intT j = 0; j < cols; j++)
        {
            if (mat[IDX(i, j, lead_dim, fmt)] != 0.) return 0;
        }
    }
    return 1;
}

int mat_cpy(DataT* in_mat, intT in_rows, intT in_cols, intT in_lead_dim, int in_fmt, DataT* out_mat, intT out_lead_dim, int out_fmt)
{
    //copy a matrix or submatrix in another matrix or submatrix of the same dimension. 

    intT in_idx,out_idx;

    //TODO add error check

    for (intT i = 0; i < in_rows; i++)
    {
        for (intT j = 0; j < in_cols; j++)
        {
            in_idx  =   IDX(i, j, in_lead_dim, in_fmt);
            out_idx =   IDX(i, j, out_lead_dim, out_fmt);

            //TODO add out of bounds check
            out_mat[out_idx] = in_mat[in_idx];
        }
    }

    return 0;
}

int random_mat(DataT* mat, intT rows, intT cols, float sparsity)
{
    //create a random matrix with given sparsity and unitary entries;

    intT nzs = (intT) (sparsity * rows * cols);
    svd entries = svd(rows * cols, 0);
    std::fill(entries.begin(), entries.begin() + nzs, 1.);

    std::random_shuffle(entries.begin(), entries.end());
    std::copy(entries.begin(), entries.end(), mat);

    return 0;

}

int equal(intT rows, intT cols, DataT* A, intT lead_A, int fmt_A, DataT* B, intT lead_B, int fmt_B, DataT eps)
{
    for (intT i = 0; i < rows; i++)
    {
        for (intT j = 0; j < cols; j++)
        {
            intT idx_A = IDX(i, j, lead_A, fmt_A);
            intT idx_B = IDX(i, j, lead_B, fmt_B);
            if (std::abs(A[idx_A] - B[idx_B]) > eps) return 0;
        }
    }
    return 1;
}

int random_sparse_blocks_mat(DataT *mat, intT rows, intT cols, int fmt, intT block_size, float block_sparsity, float block_entries_sparsity) 
{

    if ((rows % block_size != 0) or (cols % block_size != 0))
    {
        std::cout << "ERROR: matrix dimension must be a multiple of block size" << std::endl;
        return 1;
    }

    intT n_blocks = (intT)rows * cols / (block_size * block_size);     //total number of blocks
    intT nzblocks = (intT)(block_sparsity * n_blocks);              //total number of nonzero blocks

    std::fill(mat, mat + rows * cols, 0);
    svi blocks = svi(n_blocks, 0);              //will store 0 unless a block is nonzero;
    std::fill(blocks.begin(), blocks.begin() + nzblocks, 1);    //make nzblocks blocks nonzero;
    std::random_shuffle(blocks.begin(), blocks.end());          //put the nonzero blocks in random positions

    intT mat_lead_dim = (fmt == 0) ? cols : rows;

    //put nonzerovalues in the mat
    for (intT i = 0; i < rows; i += block_size) {//iterate through block rows
        intT ib = i / block_size;
        for (intT j = 0; j < cols; j += block_size) { //iterate through block columns
            intT jb = j / block_size;

            if (blocks[ib * (cols / block_size ) + jb] != 0) {
                //if block is nonempty, put random values in it;

                DataT tmp_block[block_size*block_size] = { 0 }; //temporary block
                random_mat(tmp_block, block_size, block_size, block_entries_sparsity); //make a random block with given sparsity

                intT block_idx = IDX(ib * block_size, jb * block_size, mat_lead_dim, fmt); //find starting position of the block in mat
                mat_cpy(tmp_block, block_size, block_size, block_size, 0, mat + block_idx, mat_lead_dim, fmt); //copy entries from tmp_block to mat
            }
        }
    }


}

int random_sparse_blocks_mat(VBS& vbmat, intT rows, intT cols, int blocks_fmt, int entries_fmt, intT row_block_size, intT col_block_size, float block_density, float entries_density)
{
    /*
    IN: 
        rows
        cols
        blocks_fmt : the storage format of block-rows and block-columns
        entries_fmt: the storage format INSIDE blocks
        block_density : the % of nonzero blocks
        entries_density: the % of nonzero entries in each block
        row_block_size: the height of blocks;
        col_block_size: the lenght of blocks;


    OUT: 
        vbmat : the VBS matrix, now filled.
    */

    if ((rows % row_block_size != 0) or (cols % col_block_size != 0))
    {
        //TODO exception
        std::cout << "ERROR: matrix dimension must be a multiple of block size" << std::endl;
        return 1;
    }

    intT block_rows = std::ceil((float)rows / row_block_size);
    intT block_cols = std::ceil((float)cols / col_block_size);
    intT size_of_block = row_block_size * col_block_size;
    intT n_blocks = block_rows * block_cols;
    
    intT* row_part = new intT[block_rows + 1];
    intT* col_part = new intT[block_cols + 1];
    partition(row_part, 0, rows, row_block_size); //row and block partition for creating the VBS
    partition(col_part, 0, cols, col_block_size);

    //DETERMINE THE NZ BLOCK STRUCTURE
    //(this could be probably made to use less memory by extracting indices instead of permuting the vector)
    intT nz_blocks = (intT)(block_density * block_rows * block_cols);
    svi blocks = svi(n_blocks, 0);              //will store 0 unless a block is nonzero;
    std::fill(blocks.begin(), blocks.begin() + nz_blocks, 1);    //make nzblocks blocks nonzero;
    std::random_shuffle(blocks.begin(), blocks.end());          //put the nonzero blocks in random positions

    intT main_dim = (blocks_fmt == 0)? block_rows : block_cols;
    intT compressed_dim = (blocks_fmt == 0) ? block_cols : block_rows;
    intT nz_tot = nz_blocks * size_of_block;

    init_VBS(vbmat, block_rows, row_part, block_cols, col_part, blocks_fmt, entries_fmt);

    vbmat.nztot = nz_tot;
    
    vbmat.mab = new DataT[nz_tot]{ 0 };
    vbmat.jab = new intT[nz_blocks];
    
    intT nz_in_block = std::ceil(entries_density * size_of_block);
    //saves the indices of nonzero blocks into jab; saves the number of nz blocks per row (or col) into vbmat.nzcount;
    intT b = 0;
    intT jab_count = 0;
    DataT* mab_idx = vbmat.mab; //the mab array is filled up to this pointer
    for (intT i = 0; i < main_dim; i++)
    {
        intT nzcount = 0;
        for (intT j = 0; j < compressed_dim; j++)
        {
            if (blocks[b] != 0)
            {
                std::fill(mab_idx, mab_idx + nz_in_block, 1); //fill part of the block with ones
                std::random_shuffle(mab_idx, mab_idx + size_of_block); //shuffle the block
                vbmat.jab[jab_count] = j; //save the index of the block in jab
                nzcount += 1; //keep the count of nonzero blocks
                jab_count += 1;
                mab_idx += size_of_block; // jump to next block on mab 
            }
            b++;
        }

        vbmat.nzcount[i] = nzcount;
    }

    return 0;
}

int matprint(DataT* mat, intT rows, intT cols, intT lead_dim, int fmt, bool struct_only)
{
    for (intT i = 0; i < rows; i++)
    {
        for (intT j = 0; j < cols; j++)
        {
            intT idx = IDX(i, j, lead_dim, fmt);
            std::cout << print_mat_val(mat[idx], struct_only) << " ";
        }

        std::cout << std::endl;
    }

    return 0;
}

int matprint(DataT* mat, intT rows, intT* row_part, intT row_blocks, intT cols, intT* col_part, intT col_blocks, intT lead_dim,  int fmt, bool struct_only)
{
    for (intT ib = 0; ib < row_blocks; ib++)
    {

        for (intT i = row_part[ib]; i < row_part[ib + 1]; i++)
        {
            for (intT jb = 0; jb < col_blocks; jb++)
            {
                for (intT j = col_part[jb]; j < col_part[jb + 1]; j++)
                {
                    intT idx = IDX(i, j, lead_dim, fmt);
                    std::cout << print_mat_val(mat[idx], struct_only) << " ";
                }

                std::cout << " ";

            }

            std::cout << std::endl;
        }

        std::cout << std::endl;
        
    }

    return 0;
}

int arr_print(intT* arr, intT len)
{
    for (intT i = 0; i < len; i++)
    {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;
}

int sort_permutation(intT* perm, intT* arr, intT n)
/*
returns the permutation that would sort arrary arr

    IN: perm, array of length n
        arr, array of length n
        n

   OUT: perm, the permutation that would sort arr
*/
{

    for (intT i = 0; i < n; i++) {
        perm[i] = i; //sorted indices
    }
    //sort indices in perm by hash value
    std::sort(perm, perm + n,
        [&](const intT& a, const intT& b) {return (arr[a] < arr[b]); }
    );
}

int partition(intT* arr, intT start, intT end, int step)
{
    if (step <= 0)
    {
        std::cout << "Error: step must be positive" << std::endl;
        return(0);
    }
    intT val = start;
    intT i = 0;
    while (val < end)
    {
        arr[i] = val;
        val += step;
        i++;
    }
    arr[i] = end;
}

int randperm(intT* arr, intT len)
{
    for (intT i = 0; i < len; i++)
    {
        arr[i] = i;
    }
    std::random_shuffle(arr, arr + len);
}


//TODO check this function, why return array?
intT* rand_partition(intT* part, intT len, intT blocks)
{
    intT* arr = new intT[len];
    for (intT i = 0; i < blocks; i++)
    {
        arr[i] = i;
    }
    std::random_shuffle(arr + 1, arr + len);
    std::sort(arr + 1, arr + blocks + 1);
    std::copy(arr, arr + blocks + 1, part);
    return arr;
}



//VBS utilities

int cleanVBS(VBS& vbmat)
{

    if (vbmat.nzcount) delete[] vbmat.nzcount;

    if (vbmat.jab) delete[] vbmat.jab;

    if (vbmat.mab) delete[] vbmat.mab;

    if (vbmat.row_part) delete[] vbmat.row_part;

    if (vbmat.col_part) delete[] vbmat.col_part;

    return 0;
}


//NEW
//NOT TESTED YET
int init_VBS(VBS& vbmat, intT block_rows, intT* row_part, intT block_cols, intT* col_part, int blocks_fmt, int entries_fmt)
{
    intT main_dim = (blocks_fmt == 0) ? block_rows : block_cols;
    intT compressed_dim = (blocks_fmt == 0) ? block_cols : block_rows;

    vbmat.entries_fmt = entries_fmt;
    vbmat.blocks_fmt = blocks_fmt;

    vbmat.block_rows = block_rows;
    vbmat.block_cols = block_cols;
    
    vbmat.row_part = new intT[block_rows + 1];
    std::copy(row_part, row_part + block_rows + 1, vbmat.row_part);
    
    vbmat.col_part = new intT[block_cols + 1];
    std::copy(col_part, col_part + block_cols + 1, vbmat.col_part);

    vbmat.nzcount = new intT[main_dim];


    //TODO: add conformity checks and return errors
    return 0;
}

int  convert_to_VBS(DataT* mat, intT mat_rows, intT mat_cols, int mat_fmt, VBS& vbmat, intT block_rows, intT* row_part, intT block_cols, intT *col_part, int vbmat_blocks_fmt, int vbmat_entries_fmt, int no_zero_mode)
{

    vbmat.block_rows = block_rows;
    vbmat.block_cols = block_cols;

    vbmat.blocks_fmt = vbmat_blocks_fmt;
    vbmat.entries_fmt = vbmat_entries_fmt;

    vbmat.row_part = new intT[block_rows + 1];
    vbmat.col_part = new intT[block_cols + 1];
    std::copy(row_part, row_part + block_rows + 1, vbmat.row_part);
    std::copy(col_part, col_part + block_cols + 1, vbmat.col_part);

    intT vbmat_main_dim, vbmat_compressed_dim;
    intT *b_main_ptr, *b_second_ptr;

    if (vbmat.blocks_fmt == 0)
    {
        //if Compressed Sparse Row
        vbmat_main_dim = vbmat.block_rows;
        vbmat_compressed_dim = vbmat.block_cols;

        // assign row and column pointers respectively
        // to counters for main and compressed dimension
        b_main_ptr = row_part;
        b_second_ptr = col_part;

    }
    else
    {
        //if Compressed Sparse Columns
        vbmat_main_dim = vbmat.block_cols;
        vbmat_compressed_dim = vbmat.block_rows;

        // assign column and row pointers respectively
        // to counters for main and compressed dimension
        b_main_ptr = col_part;
        b_second_ptr = row_part;
    }

    intT main_pos, main_block_dim, second_pos, second_block_dim;
    intT row, col, row_block_dim, col_block_dim;
    intT total_nonzero_entries = 0; 

    vbmat.nzcount = new intT[vbmat_main_dim]; //initialize nzcount, which stores number of nonzero blocks for each block-row (-column)
    intT mat_leading_dim = mat_fmt == 0 ? mat_cols : mat_rows;

    intT matpos;
    svi jab;
    svd mab;


    //FIND BLOCK STRUCTURE--------------------------------------------------------------
    for (intT i = 0; i < vbmat_main_dim; i++)        //loops trough main block dimension
    {
        main_pos = b_main_ptr[i];
        main_block_dim = b_main_ptr[i + 1] - main_pos;
        vbmat.nzcount[i] = 0;

        for (intT j = 0; j < vbmat_compressed_dim; j++)     //loops through compressed block dimension
        {

            second_pos = b_second_ptr[j];
            second_block_dim = b_second_ptr[j + 1] - second_pos;

            if (vbmat.blocks_fmt == 0)
            {
                row = main_pos;
                col = second_pos;
                row_block_dim = main_block_dim;
                col_block_dim = second_block_dim;
            }
            else
            {
                row = second_pos;
                col = main_pos;
                row_block_dim = second_block_dim;
                col_block_dim = main_block_dim;
            }

            matpos = IDX(row, col, mat_leading_dim, mat_fmt);    //find starting index of block in matrix
            if (!is_empty(mat + matpos, row_block_dim, col_block_dim, mat_leading_dim, mat_fmt) or no_zero_mode == 1)        //check if block is non-empty
            {
                vbmat.nzcount[i] += 1;  //one more nonzero block on the compressed dimension
                total_nonzero_entries += main_block_dim * second_block_dim;
                jab.push_back(j);       //store index of nonzero block
            }
        }
    }
    //----------------------------------------------------------------------------------

    vbmat.nztot = total_nonzero_entries;
    vbmat.jab = new intT[jab.size()];
    vbmat.mab = new DataT[total_nonzero_entries];

    std:copy(jab.begin(), jab.end(), vbmat.jab);

    intT mat_idx = 0; //keeps reading position for mat
    intT vbmat_idx = 0; //keeps writing position for vbmat 
    intT jab_count = 0;

    //COPY VALUES from mat to vbmat ------------------------------------------------------
    for (intT i = 0; i < vbmat_main_dim; i++)
    {
        main_pos = b_main_ptr[i];
        main_block_dim = b_main_ptr[i + 1] - main_pos;

     

        for (intT nzs = 0; nzs < vbmat.nzcount[i]; nzs++)
        {
            intT j = jab[jab_count];

            second_pos = b_second_ptr[j];
            second_block_dim = b_second_ptr[j + 1] - second_pos;

            if (vbmat.blocks_fmt == 0)
            {
                row = main_pos;
                col = second_pos;
                row_block_dim = main_block_dim;
                col_block_dim = second_block_dim;
            }
            else
            {
                row = second_pos;
                col = main_pos;
                row_block_dim = second_block_dim;
                col_block_dim = main_block_dim;
            }

            mat_idx = IDX(row, col, mat_leading_dim, mat_fmt); //find starting index of block in matrix

            intT block_leading_dim = (vbmat.entries_fmt == 0) ? col_block_dim : row_block_dim;

            mat_cpy(mat + mat_idx, row_block_dim, col_block_dim, mat_leading_dim, mat_fmt, vbmat.mab + vbmat_idx, block_leading_dim, vbmat_entries_fmt); //write block from mat to vbmat.mab
            
            vbmat_idx += main_block_dim * second_block_dim;
            jab_count++;
        }
    }

    return 0;
}

int convert_to_mat(const VBS& vbmat, DataT* out_mat, int out_mat_fmt)
{
    //input:
    //  vmat: a VBS matrix
    //  out_mat: an array of the proper dimension, filled with 0s; 

    //determine out_mat dimensions-------------------------

    intT out_mat_rows = vbmat.row_part[vbmat.block_rows];
    intT out_mat_cols = vbmat.col_part[vbmat.block_cols];
    intT mat_leading_dim = out_mat_fmt == 0 ? out_mat_cols : out_mat_rows;
    //-----------------------------------------------------

    intT main_pos, main_block_dim, second_pos, second_block_dim;
    intT row, col, row_block_dim, col_block_dim;

    intT mat_idx = 0; //keeps writing position for mat
    intT vbmat_idx = 0; //keeps reading position for vbmat 

    intT vbmat_main_dim = vbmat.main_dim();
    intT vbmat_compressed_dim = vbmat.compressed_dim();
    intT *b_main_ptr = vbmat.main_ptr();
    intT *b_second_ptr = vbmat.second_ptr();

    intT* jab = vbmat.jab;

    intT nz_tot = 0; //counter for total nonzero blocks
    //COPY VALUES from vbmat to mat ------------------------------------------------------
    for (intT i = 0; i < vbmat_main_dim; i++)
    {
        main_pos = b_main_ptr[i];
        main_block_dim = b_main_ptr[i + 1] - main_pos;

        for (intT nzs = 0; nzs < vbmat.nzcount[i]; nzs++) //iterate for all nonzero block in row i
        {
            intT j = jab[nz_tot]; //column of current non-zero block
            nz_tot += 1; //count one nonzero block

            second_pos = b_second_ptr[j];
            second_block_dim = b_second_ptr[j + 1] - second_pos;

            if (vbmat.blocks_fmt == 0)
            {
                row = main_pos;
                col = second_pos;
                row_block_dim = main_block_dim;
                col_block_dim = second_block_dim;
            }
            else
            {
                row = second_pos;
                col = main_pos;
                row_block_dim = second_block_dim;
                col_block_dim = main_block_dim;
            }

            mat_idx = IDX(row, col, mat_leading_dim, out_mat_fmt); //find starting index of block in matrix

            intT block_leading_dim = (vbmat.entries_fmt == 0) ? col_block_dim : row_block_dim;

            mat_cpy(vbmat.mab + vbmat_idx, row_block_dim, col_block_dim, block_leading_dim, vbmat.entries_fmt, out_mat + mat_idx, mat_leading_dim, out_mat_fmt); //write block from vbmat.mab to mat

            vbmat_idx += row_block_dim * col_block_dim; //update vbmat.mab reading index
        }
    }
    //------------------------------------------------------------------------------------

    return 0;

}

int convert_to_VBS(const CSR& cmat, VBS& vbmat, intT block_rows, intT* row_part, intT block_cols, intT* col_part, int vbmat_block_fmt, int vbmat_entries_fmt)
{

    intT mat_rows = cmat.rows;
    intT mat_cols = cmat.cols;
    
    intT cmat_main_dim = cmat.fmt == 0 ? cmat.rows : cmat.cols;
    intT cmat_second_dim = cmat.fmt == 0 ? cmat.cols : cmat.rows;

    init_VBS(vbmat, block_rows, row_part, block_cols, col_part, vbmat_block_fmt, vbmat_entries_fmt);

    intT vbmat_main_dim = vbmat_block_fmt == 0 ? block_rows : block_cols;
    intT vbmat_second_dim = vbmat_block_fmt == 0 ? block_cols : block_rows;

    intT total_blocks = block_rows * block_cols;
    
    intT** blocks_bookmark = new intT*[vbmat_main_dim]; //will identify nonzero blocks;
    for (intT i = 0; i < vbmat_main_dim; i++)
    {
        blocks_bookmark[i] = new intT[vbmat_second_dim];
        for (intT j = 0; j < vbmat_second_dim; j++)
        {
            blocks_bookmark[i][j] = -1;
        }
    }

    //pointers for proper iteration depending on fmt
    intT current_block_row = 0, current_block_col = 0;
    intT* current_main_pos = new intT();
    intT* current_second_pos = new intT();
    if (vbmat_block_fmt == 0)
    {
        current_main_pos = &current_block_row;
        current_second_pos = &current_block_col;
    }
    else
    {
        current_main_pos = &current_block_col;
        current_second_pos = &current_block_row;
    }


    //bookmarks nonzero blocks
    for (intT i = 0; i < cmat_main_dim; i++)
    {
        current_block_row = 0;
        current_block_col = 0;

        for (intT nzs = 0; nzs < cmat.nzcount[i]; nzs++)
        {

            intT j = cmat.ja[i][nzs];
            //find vbmat block;
            intT row = cmat.fmt == 0 ? i : j;
            intT col = cmat.fmt == 0 ? j : i;


            while (row >= row_part[current_block_row])
            {
                current_block_row++;
            }

            while (col >= col_part[current_block_col])
            {
                current_block_col++;
            }

            current_block_row--;
            current_block_col--;

            //flag the bookmark position (nonzero block)
            blocks_bookmark[*current_main_pos][*current_second_pos] = -2;
        }

    }
    //DEBUG 
    int test = 0;


    //counts total nonzero area and total nonzero blocks
    intT total_area = 0;
    intT total_nz_blocks = 0;
    for (intT ib = 0; ib < vbmat_main_dim; ib++)
    {
        vbmat.nzcount[ib] = 0;
        for (intT jb = 0; jb < vbmat_second_dim; jb++)
        {
            intT block_row = vbmat.blocks_fmt == 0 ? ib : jb;
            intT block_col = vbmat.blocks_fmt == 0 ? jb : ib;
            if (blocks_bookmark[ib][jb] == -2)
            {
                blocks_bookmark[ib][jb] = total_area;
                total_area += (row_part[block_row + 1] - row_part[block_row]) * (col_part[block_col + 1] - col_part[block_col]);
                vbmat.nzcount[ib] += 1;
                total_nz_blocks += 1;
            }
        }
    }

    //assign the values to the matrix
    vbmat.mab = new DataT[total_area]{ 0 };
    vbmat.nztot = total_area;
    vbmat.jab = new intT[total_nz_blocks];
   
    

    //fill jab (nonzero blocks positions)
    total_nz_blocks = 0;
    for (intT ib = 0; ib < vbmat_main_dim; ib++)
    {
        for (intT jb = 0; jb < vbmat_second_dim; jb++)
        {
            if (blocks_bookmark[ib][jb] != -1)
            {
                vbmat.jab[total_nz_blocks] = jb;
                total_nz_blocks += 1;
            }
        }
    }


    //fill vbmat nonzeros with entries from cmat
    for (intT i = 0; i < cmat_main_dim; i++)
    {
        current_block_row = 0;
        current_block_col = 0;

        for (intT nzs = 0; nzs < cmat.nzcount[i]; nzs++)
        {
            intT j = cmat.ja[i][nzs];
            //find vbmat block;
            intT row = cmat.fmt == 0 ? i : j;
            intT col = cmat.fmt == 0 ? j : i;

            while (row >= row_part[current_block_row])
            {
                current_block_row++;
            }

            while (col >= col_part[current_block_col])
            {
                current_block_col++;
            }
            current_block_row--;
            current_block_col--;


            //put the value in the correct block at the correct position

            intT block_start = blocks_bookmark[*current_main_pos][*current_second_pos];
            intT block_rows = row_part[current_block_row + 1] - row_part[current_block_row];
            intT block_cols = col_part[current_block_col + 1] - col_part[current_block_col];
            intT relative_row_pos = row - row_part[current_block_row];
            intT relative_col_pos = col - col_part[current_block_col];

            intT block_leading_dim = vbmat.entries_fmt == 0 ? block_cols : block_rows;
            intT block_idx = block_start + IDX(relative_row_pos, relative_col_pos, block_leading_dim, vbmat.entries_fmt);
            vbmat.mab[block_idx] = cmat.ma[i][nzs];
        }

    }

    for (intT i = 0; i < vbmat_main_dim; i++)
    {
        delete[] blocks_bookmark[i];
    }
    delete[] blocks_bookmark;
    return 0;
}

int convert_to_VBS(const CSR& cmat, VBS& vbmat, intT row_block_size, intT col_block_size, int vbmat_block_fmt, int vbmat_entries_fmt)
{

    intT mat_rows = cmat.rows;
    intT mat_cols = cmat.cols;

    intT cmat_main_dim = cmat.fmt == 0 ? cmat.rows : cmat.cols;
    intT cmat_second_dim = cmat.fmt == 0 ? cmat.cols : cmat.rows;


    intT block_rows = mat_rows / row_block_size;
    intT block_cols = mat_cols / col_block_size;
    
    intT total_blocks = block_rows * block_cols;

    intT* row_part = new intT[block_rows + 1]; //partitions have one element more for the rightmost border.
    intT* col_part = new intT[block_cols + 1];
    partition(row_part, 0, mat_rows, row_block_size); //row and column partitions (TODO make it work when block_size does not divide rows)
    partition(col_part, 0, mat_cols, col_block_size);

    init_VBS(vbmat, block_rows, row_part, block_cols, col_part, vbmat_block_fmt, vbmat_entries_fmt);

    intT vbmat_main_dim = vbmat_block_fmt == 0 ? block_rows : block_cols;
    intT vbmat_second_dim = vbmat_block_fmt == 0 ? block_cols : block_rows;

    intT** blocks_bookmark = new intT * [vbmat_main_dim]; //will identify nonzero blocks;
    for (intT i = 0; i < vbmat_main_dim; i++)
    {
        blocks_bookmark[i] = new intT[vbmat_second_dim];
        for (intT j = 0; j < vbmat_second_dim; j++)
        {
            blocks_bookmark[i][j] = -1;
        }
    }

    //pointers for proper iteration depending on fmt
    intT current_block_row = 0, current_block_col = 0;
    intT* current_main_pos = new intT();
    intT* current_second_pos = new intT();
    if (vbmat_block_fmt == 0)
    {
        current_main_pos = &current_block_row;
        current_second_pos = &current_block_col;
    }
    else
    {
        current_main_pos = &current_block_col;
        current_second_pos = &current_block_row;
    }


    //bookmarks nonzero blocks
    for (intT i = 0; i < cmat_main_dim; i++)
    {
        current_block_row = 0;
        current_block_col = 0;

        for (intT nzs = 0; nzs < cmat.nzcount[i]; nzs++)
        {

            intT j = cmat.ja[i][nzs];
            //find vbmat block;
            intT row = cmat.fmt == 0 ? i : j;
            intT col = cmat.fmt == 0 ? j : i;


            while (row >= row_part[current_block_row])
            {
                current_block_row++;
            }

            while (col >= col_part[current_block_col])
            {
                current_block_col++;
            }

            current_block_row--;
            current_block_col--;

            //flag the bookmark position (nonzero block)
            blocks_bookmark[*current_main_pos][*current_second_pos] = -2;
        }

    }
    //DEBUG 
    int test = 0;


    //counts total nonzero area and total nonzero blocks
    intT total_area = 0;
    intT total_nz_blocks = 0;
    for (intT ib = 0; ib < vbmat_main_dim; ib++)
    {
        vbmat.nzcount[ib] = 0;
        for (intT jb = 0; jb < vbmat_second_dim; jb++)
        {
            intT block_row = vbmat.blocks_fmt == 0 ? ib : jb;
            intT block_col = vbmat.blocks_fmt == 0 ? jb : ib;
            if (blocks_bookmark[ib][jb] == -2)
            {
                blocks_bookmark[ib][jb] = total_area;
                total_area += (row_part[block_row + 1] - row_part[block_row]) * (col_part[block_col + 1] - col_part[block_col]);
                vbmat.nzcount[ib] += 1;
                total_nz_blocks += 1;
            }
        }
    }

    //assign the values to the matrix
    vbmat.mab = new DataT[total_area]{ 0 };
    vbmat.nztot = total_area;
    vbmat.jab = new intT[total_nz_blocks];



    //fill jab (nonzero blocks positions)
    total_nz_blocks = 0;
    for (intT ib = 0; ib < vbmat_main_dim; ib++)
    {
        for (intT jb = 0; jb < vbmat_second_dim; jb++)
        {
            if (blocks_bookmark[ib][jb] != -1)
            {
                vbmat.jab[total_nz_blocks] = jb;
                total_nz_blocks += 1;
            }
        }
    }


    //fill vbmat nonzeros with entries from cmat
    for (intT i = 0; i < cmat_main_dim; i++)
    {
        current_block_row = 0;
        current_block_col = 0;

        for (intT nzs = 0; nzs < cmat.nzcount[i]; nzs++)
        {
            intT j = cmat.ja[i][nzs];
            //find vbmat block;
            intT row = cmat.fmt == 0 ? i : j;
            intT col = cmat.fmt == 0 ? j : i;

            while (row >= row_part[current_block_row])
            {
                current_block_row++;
            }

            while (col >= col_part[current_block_col])
            {
                current_block_col++;
            }
            current_block_row--;
            current_block_col--;


            //put the value in the correct block at the correct position

            intT block_start = blocks_bookmark[*current_main_pos][*current_second_pos];
            intT block_rows = row_part[current_block_row + 1] - row_part[current_block_row];
            intT block_cols = col_part[current_block_col + 1] - col_part[current_block_col];
            intT relative_row_pos = row - row_part[current_block_row];
            intT relative_col_pos = col - col_part[current_block_col];

            intT block_leading_dim = vbmat.entries_fmt == 0 ? block_cols : block_rows;
            intT block_idx = block_start + IDX(relative_row_pos, relative_col_pos, block_leading_dim, vbmat.entries_fmt);
            vbmat.mab[block_idx] = cmat.ma[i][nzs];
        }

    }

    for (intT i = 0; i < vbmat_main_dim; i++)
    {
        delete[] blocks_bookmark[i];
    }
    delete[] blocks_bookmark;
    return 0;
}


int matprint(const VBS& vbmat)
{

    intT rows = vbmat.row_part[vbmat.block_rows];    
    intT cols = vbmat.col_part[vbmat.block_cols];

    DataT* temp_mat = new DataT[rows * cols]{0};

    convert_to_mat(vbmat, temp_mat, 0);
    
    matprint(temp_mat, rows, vbmat.row_part, vbmat.block_rows, cols, vbmat.col_part, vbmat. block_cols, cols, 0);
    delete[] temp_mat;

    return 0;
}


//CSR utilities

int cleanCSR(CSR& cmat)
{
    /*----------------------------------------------------------------------
    | Free up memory allocated for CSR structs.
    |----------------------------------------------------------------------
    | on entry:
    |==========
    |  cmat  =  a CSR struct.
    |--------------------------------------------------------------------*/
    /*   */

    if (cmat.rows + cmat.cols <= 1) return 0;

    intT main_dim;
    main_dim = cmat.fmt == 0 ? cmat.rows : cmat.cols;

    for (intT i = 0; i < main_dim; i++) {
        if (cmat.nzcount[i] > 0) {
            if (cmat.ma) delete[] cmat.ma[i];
            delete[] cmat.ja[i];
        }
    }
    if (cmat.ma) delete[] cmat.ma;
    delete[] cmat.ja;
    delete[] cmat.nzcount;

    return 0;
}

int convert_to_mat(const CSR& cmat, DataT* out_mat, int out_mat_fmt)
{
    //input:
    //  vmat: a CSR matrix
    //  out_mat: an array of the proper dimension, filled with 0s; 

    intT main_dim;
    intT compressed_dim;
    intT i; //counter for main dimension
    intT j; //counter for compressed dimension

    intT* row_ptr; //pointer to current block_row
    intT* col_ptr; //pointer to current block_column

    if (cmat.fmt == 0)
    {
        //if Compressed Sparse Row
        main_dim = cmat.rows;
        compressed_dim = cmat.cols;

        // assign row and column pointers respectively
        // to counters for main and compressed dimension
        row_ptr = &i;
        col_ptr = &j;

    }
    else
    {
        //if Compressed Sparse Columns
        main_dim = cmat.cols;
        compressed_dim = cmat.rows;

        // assign column and row pointers respectively
        // to counters for main and compressed dimension
        row_ptr = &j;
        col_ptr = &i;
    }

    intT out_lead_dim = out_mat_fmt == 0 ? cmat.cols : cmat.rows;
    
    //FILL THE OUTPUT MATRIX
    for (i = 0; i < main_dim; i++)
    {
        for (intT nzs = 0; nzs < cmat.nzcount[i]; nzs++) 
        {
            j = cmat.ja[i][nzs]; //find column (row) index of next nonzero element
            DataT elem = cmat.ma[i][nzs]; //value of that element;

            intT mat_idx = IDX(*row_ptr, *col_ptr, out_lead_dim, out_mat_fmt);
            out_mat[mat_idx] = elem;
        }
    }

    return 0;

}

int convert_to_CSR(const DataT* in_mat, intT mat_rows, intT mat_cols, int mat_fmt, CSR& cmat, int cmat_fmt)
{

    cmat.rows = mat_rows;
    cmat.cols = mat_cols;
    cmat.fmt = cmat_fmt;

    intT cmat_main_dim;
    intT cmat_compressed_dim;
    intT i; //counter for main dimension
    intT j; //counter for compressed dimension

    intT* cmat_row_ptr; //pointer to current block_row
    intT* cmat_col_ptr; //pointer to current block_column


    if (cmat_fmt == 0)
    {
        //if Compressed Sparse Row
        cmat_main_dim = cmat.rows;
        cmat_compressed_dim = cmat.cols;

        // assign row and column pointers respectively
        // to counters for main and compressed dimension
        cmat_row_ptr = &i;
        cmat_col_ptr = &j;

    }
    else
    {
        //if Compressed Sparse Columns
        cmat_main_dim = cmat.cols;
        cmat_compressed_dim = cmat.rows;

        // assign column and row pointers respectively
        // to counters for main and compressed dimension
        cmat_row_ptr = &j;
        cmat_col_ptr = &i;
    }

    cmat.nzcount = new intT[cmat_main_dim];
    cmat.ja = new intT* [cmat_main_dim];
    cmat.ma = new DataT* [cmat_main_dim];

    intT* tmp_ja = new intT[cmat_compressed_dim];
    DataT* tmp_ma = new DataT[cmat_compressed_dim];

    intT mat_lead_dim = mat_fmt == 0 ? cmat.cols : cmat.rows;

    for (i = 0; i < cmat_main_dim; i++) //loop through main dimension
    {
        intT nzs = 0; //non-zero entries in this row (column) 
        for (j = 0; j < cmat_compressed_dim; j++) //scan compressed dimension for nonzero entries
        {
            intT mat_idx = IDX(*cmat_row_ptr, *cmat_col_ptr, mat_lead_dim, mat_fmt);
            if (in_mat[mat_idx] != 0) //if the entry is nonzero add its idx to ja and its value to ma
            {
                tmp_ja[nzs] = j;
                tmp_ma[nzs] = in_mat[mat_idx];
                nzs++;
            }
        }

        cmat.nzcount[i] = nzs; 
        cmat.ja[i] = new intT[nzs];
        cmat.ma[i] = new DataT[nzs];

        std::copy(tmp_ja, tmp_ja + nzs, cmat.ja[i]); 
        std::copy(tmp_ma, tmp_ma + nzs, cmat.ma[i]);

    }

    delete[] tmp_ma;
    delete[] tmp_ja;

    return 0;
}


int convert_to_CSR(const VBS& vbmat, CSR& cmat, int csr_fmt = 0)
{

    //ONLY ACCEPTS fmt = 0 for the moment;
    cmat.rows = vbmat.rows();
    cmat.cols = vbmat.cols();
    cmat.fmt = csr_fmt;

    cmat.nzcount = new intT[cmat.rows]{ 0 };
    cmat.ja = new intT * [cmat.rows];
    cmat.ma = new DataT * [cmat.rows];

    intT main_pos;
    intT main_block_dim;
    intT second_block_dim;
    intT row_start; //first row of a block
    intT col_start; //first col of a block
    intT jb; 
    intT* jab = vbmat.jab;
    DataT* mab = vbmat.mab;

    for (intT ib = 0; ib < vbmat.main_dim(); ib++) // loop through vbmat main dim
    {

        main_block_dim = vbmat.blocks_fmt ? vbmat.block_width(ib) : vbmat.block_height(ib);

        for (intT nzs = 0; nzs < vbmat.nzcount[ib]; nzs++) //loop through nonzero blocks on this row
        {
            jb = jab[0];
            jab++;
            second_block_dim = vbmat.blocks_fmt ? vbmat.block_height(jb) : vbmat.block_width(jb);

            row_start = vbmat.blocks_fmt ? vbmat.row_part[jb] : vbmat.row_part[ib];
            col_start = vbmat.blocks_fmt ? vbmat.col_part[ib] : vbmat.col_part[jb];

            intT row;
            intT col;
            for (intT i = 0; i < main_block_dim; i++) //loop through entries in the block
            {
                for (intT j = 0; j < second_block_dim; j++)
                {
                    if (mab[0] != 0)
                    {
                        row = vbmat.entries_fmt == 0 ? i + row_start : j + row_start;
                        col = vbmat.entries_fmt == 0 ? j + col_start : i + col_start;
                        cmat.nzcount[row]++; //found a nonzero in this row. Added to the count
                    }
                    mab++; //procede forward in the vbmat entries array
                }
            }
        }
    }

    for (intT i = 0; i < cmat.rows; i++)
    {
        if (cmat.nzcount[i])
        {
            cmat.ja[i] = new intT[cmat.nzcount[i]];
            cmat.ma[i] = new DataT[cmat.nzcount[i]];
        }
    }

    intT* current_nz_count = new intT[cmat.rows]{ 0 };
    
    jab = vbmat.jab; //reset pointers
    mab = vbmat.mab;
    for (intT ib = 0; ib < vbmat.main_dim(); ib++) // loop through vbmat main dim
    {

        main_block_dim = vbmat.blocks_fmt ? vbmat.block_width(ib) : vbmat.block_height(ib);

        for (intT nzs = 0; nzs < vbmat.nzcount[ib]; nzs++) //loop through nonzero blocks on this row
        {
            jb = jab[0];
            jab++;

            second_block_dim = vbmat.blocks_fmt ? vbmat.block_height(jb) : vbmat.block_width(jb);

            row_start = vbmat.blocks_fmt ? vbmat.row_part[jb] : vbmat.row_part[ib];
            col_start = vbmat.blocks_fmt ? vbmat.col_part[ib] : vbmat.col_part[jb];

            intT row;
            intT col;
            for (intT i = 0; i < main_block_dim; i++) //loop through entries in the block
            {
                for (intT j = 0; j < second_block_dim; j++)
                {
                    if (mab[0] != 0)
                    {
                        row = vbmat.entries_fmt == 0 ? i + row_start : j + row_start;
                        col = vbmat.entries_fmt == 0 ? j + col_start : i + col_start;
                        cmat.ja[row][current_nz_count[row]] = col; //save the position of the nz in the array ja[row], at the last unoccupied index
                        cmat.ma[row][current_nz_count[row]] = mab[0]; //save the value of the nz in the array ma[row], at the last unoccupied index
                        current_nz_count[row]++;
                    }
                    mab++;
                }
            }
        }
    }

    return 0;
}

int convert_to_CSR_ineff(const VBS& vbmat, CSR& cmat, int csr_fmt)
{
    //TODO: if necessary, make conversion efficient;

    intT mat_rows = vbmat.row_part[vbmat.block_cols];
    intT mat_cols = vbmat.col_part[vbmat.block_cols];
    intT mat_size = mat_rows * mat_cols;
    int mat_fmt = 0;

    DataT* mat = new DataT[mat_size]{ 0 };
    convert_to_mat(vbmat, mat, mat_fmt);
    convert_to_CSR(mat, mat_rows, mat_cols, mat_fmt, cmat, csr_fmt);
    delete[] mat;

    return 0;
}

int matprint(const CSR& cmat)
{
    DataT* tmp_mat = new DataT[cmat.rows * cmat.cols]{ 0 };
    convert_to_mat(cmat, tmp_mat, 0);
    matprint(tmp_mat, cmat.rows, cmat.cols, cmat.cols, 0);
    delete[] tmp_mat;

    return 0;
}

int copy(const CSR& in_cmat, CSR& out_cmat) 
{

    out_cmat.fmt = in_cmat.fmt;
    out_cmat.rows = in_cmat.rows;
    out_cmat.cols = in_cmat.cols;

    intT main_dim = (in_cmat.fmt == 0) ? in_cmat.rows : in_cmat.cols;

    out_cmat.nzcount = new intT[main_dim];
    out_cmat.ja = new intT* [main_dim];
    out_cmat.ma = new DataT * [main_dim];
    std::copy(in_cmat.nzcount, in_cmat.nzcount + main_dim, out_cmat.nzcount);

    for (intT i = 0; i < main_dim; i++)
    {
        intT nzs = out_cmat.nzcount[i];
        out_cmat.ja[i] = new intT[nzs];
        out_cmat.ma[i] = new DataT[nzs];

        std::copy((in_cmat.ja)[i], (in_cmat.ja)[i] + nzs, (out_cmat.ja)[i]);
        std::copy((in_cmat.ma)[i], (in_cmat.ma)[i] + nzs, (out_cmat.ma)[i]);

    }

    return 0;
}

int transpose(const CSR& in_cmat, CSR& out_cmat, int new_fmt)
{
    intT in_main_dim = in_cmat.fmt == 0 ? in_cmat.rows : in_cmat.cols; //main dimension of in_mat; secondary of out_mat; 
    intT in_second_dim = in_cmat.fmt == 0 ? in_cmat.cols : in_cmat.cols; //secondary dimension of in_mat; main for out_mat;

    if (new_fmt != in_cmat.fmt)
    {

        copy(in_cmat, out_cmat);
        out_cmat.fmt == new_fmt;

    }
    else
    {
        out_cmat.fmt = in_cmat.fmt;
        out_cmat.rows = in_cmat.cols;
        out_cmat.cols = in_cmat.rows;
    }


    out_cmat.nzcount = new intT[in_second_dim];
    out_cmat.ja = new intT* [in_second_dim];
    out_cmat.ma = new DataT *[in_second_dim];
    

    //find number of nonzero elements in each secondary row (which will be main row for the transpose); 
    for (intT i = 0; i < in_main_dim; i++)
    {
        for (intT nzs = 0; nzs < in_cmat.nzcount[i]; nzs++)
        {
            intT j = in_cmat.ja[i][nzs]; //find column (row) index of next nonzero element
            out_cmat.nzcount[j] ++;
        }
    }

    intT* counter = new intT[in_second_dim]{ 0 };
   
    //initialize arrays in out_cmat
    for (intT j = 0; j < in_second_dim; j++)
    {
        out_cmat.ja[j] = new intT[out_cmat.nzcount[j]];
        out_cmat.ma[j] = new DataT[out_cmat.nzcount[j]];
    }


    intT c = 0;
    for (intT i = 0; i < in_main_dim; i++)
    {
        for (intT nzs = 0; nzs < in_cmat.nzcount[i]; nzs++)
        {
            intT j = in_cmat.ja[i][nzs]; //find in_cmat main dim index of next nonzero element
            DataT elem = in_cmat.ma[i][nzs]; //value of that element;

            c = counter[j]; //progressively fill out_cmat main dim
            out_cmat.ja[j][c] = i;
            out_cmat.ma[j][c] = elem;
        }
    }
    
    delete[] counter;

    return 0;

}

int permute_CSR(CSR& cmat, intT* perm, int dim) {
    //permutes rows (dim == 0), cols (dim == 1) or both (dim==2) of a matrix in CSR form;
    //dim = cmat.fmt will permute the main dimension;


    intT main_dim = (cmat.fmt == 0) ? cmat.rows : cmat.cols;
    intT second_dim = (cmat.fmt == 0) ? cmat.cols : cmat.rows;


    bool permute_main = (cmat.fmt == 0) ? (dim == 0) : (dim == 1);  //permute main dimension?
    bool permute_second = !permute_main;                        //permute secondary dimension?

    if (dim == 2)
    {
        if (cmat.rows != cmat.cols)
        {
            std::cout << "ERROR: matrix must be square to apply same permutation to row and column" << std::endl;
            return 1;
        }
        permute_main = true;
        permute_second = true;
    }

    if (permute_main)
    {
        int err_check = 0;

        err_check += permute(cmat.nzcount, perm, main_dim);

        err_check += permute(cmat.ja, perm, main_dim);

        err_check += permute(cmat.ma, perm, main_dim);

        if (err_check != 0)
        {
            std::cout << "an error occurred permuting the matrix" << std::endl;
            return 1;
        }
    }

    if (permute_second)
    {

        intT* idx_perm = new intT[second_dim] {0};

        for (intT j = 0; j < second_dim; j++)
        {
            idx_perm[perm[j]] = j;   //this array stores the new names of column idxs after the permutation (i.e. for permutation 3 1 0 2, it stores 2 1 3 0) 
        }

        intT* ja;
        for (intT i = 0; i < main_dim; i++)
        {
            ja = cmat.ja[i];
            intT ja_len = cmat.nzcount[i];
            if (ja_len > 0)
            {

                //change column indices to new values (given by idx_perm)
                for (intT j = 0; j < ja_len; j++)
                {
                    ja[j] = idx_perm[ja[j]]; //assign the new names to indices
                }

                intT* tmp_perm = new intT[ja_len] {0};
                sort_permutation(tmp_perm, ja, ja_len); //find correct reorder of column indices.

                permute(cmat.ja[i], tmp_perm, ja_len); //sort ja[i]

                permute(cmat.ma[i], tmp_perm, ja_len); //permute ma[i] with the same permutation;
                delete[] tmp_perm;
            }
        }

        delete[] idx_perm;
    }

    return 0;
}

intT count_nnz(CSR& cmat)
{
    intT nnz = 0;
    intT main = (cmat.fmt == 0) ? cmat.rows : cmat.cols;
    for (intT i = 0; i < main; i++)
    {
        nnz += cmat.nzcount[i];
    }
    return nnz;
}

intT count_nnz_blocks(VBS& vbmat)
{
    intT main_block_dim = (vbmat.blocks_fmt == 0) ? vbmat.block_rows : vbmat.block_cols;
    intT nz_blocks = 0;
    for (intT ib = 0; ib < main_block_dim; ib++)
    {
        nz_blocks += vbmat.nzcount[ib];
    }
    return nz_blocks;
}

//TODO get_pattern (and other pattern related functions) can be made more efficient by storing compressed blocks instead of an array of blocks

intT scalar_product(intT* pat_0, intT len_0, intT* pat_1)
{
    intT scalar = 0;

    for (intT i = 0; i < len_0; i++)
    {
        scalar += pat_0[i] * pat_1[i];
    }

    return scalar;

}

intT norm2(intT* arr, intT len)
{
    intT norm = 0;
    for (intT i = 0; i < len; i++)
    {
        norm += arr[i] * arr[i];
    }
    return norm;
}

void read_snap_format(GraphMap& gmap, std::string filename, std::string delimiter = " ")
{
    /*		Read from edgelist to a graphmap
     *		TODO Error handling
     *----------------------------------------------------------------------------
     * on entry:
     * =========
     * filename: the file were the edgelist is stored
     *----------------------------------------------------------------------------
     * on return:
     * ==========
     * gmap   = the edgelist now into GraphMap format
     *----------------------------------------------------------------------------
     */

    gmap.clear();

    std::ifstream infile;

    //TODO handle reading error
    infile.open(filename);
    std::string temp = "-1";
    intT current_node = -1, child;
    std::set<intT> emptyset;

    // Ignore comments headers
    while (infile.peek() == '#' or infile.peek() == '%') infile.ignore(2048, '\n');

    //TODO: check import success
    while (getline(infile, temp)) {

        int del_pos = temp.find(delimiter);
        int del_size = delimiter.length();

        std::string first_node_string = temp.substr(0, del_pos); //retrieve the part of the string before the delimiter
        current_node = stoi(first_node_string);
        temp.erase(0, del_pos + del_size);

        if (gmap.count(current_node) == 0) { //new source node encountered
            gmap[current_node] = emptyset;
        }

        del_pos = temp.find(delimiter);
        std::string second_node_string = temp.substr(0, del_pos); //retrieve the part of the string after the delimiter
        child = stoi(second_node_string);
        gmap[current_node].insert(child);
    }
}
//check if a graphmap is well-numbered, complete and (optional) symmetric
int isProper(const GraphMap& gmap, bool mirroring) 
{
    //returns:
    //0 : everything ok
    //1 : there is a jump in node numeration
    //2 : inexistent children
    //3 : asymmetric link  (only if mirroring is on)

    int check = 0;
    intT last_node = -1;

    for (auto const& x : gmap)
    {
        //check if there is a jump
        if (x.first > last_node + 1) { return 1; }
        last_node = x.first;

        //check node children
        for (auto const& elem : x.second)
        {
            //check if the children exists in the graph
            auto child = gmap.find(elem); //search result
            if (child == gmap.end()) {
                //inexistent children
                return 2;
            }
            //check for asymmetric links (only if mirroring is on)
            else if (mirroring) {
                if ((child->second).find(x.first) == (child->second).end()) {
                    return 3;
                }
            }
        }
    }
    return 0;
}

void MakeUndirected(GraphMap& gmap) {
    //Make a graph undirected

    for (auto x : gmap)
    {
        //check node children
        for (auto child : x.second) {
            (gmap[child]).insert(x.first);
        }
    }
}

void MakeProper(GraphMap& gmap) {
    //rename graph nodes so that they are consecutive integers 
    //Quite costly. Probably worth checking with isProper first
    std::map<intT, intT> new_name;
    intT count = -1;
    intT name = -1;
    std::set<intT> tempset;

    //find correct names
    for (auto parent : gmap)
    {
        count++;
        name = parent.first;
        new_name[name] = count;
    }

    //change target names
    for (auto parent : gmap)
    {

        tempset.clear();
        for (auto child : parent.second) {
            tempset.insert(new_name[child]);
        }
        gmap[parent.first] = tempset;
    }

    //change source names
    for (auto n : new_name)
    {
        if (n.first != n.second) {
            gmap[n.second] = gmap[n.first];
            gmap.erase(n.first);
        }
    }
}

void write_snap_format(GraphMap& gmap, std::string filename) {
    std::ofstream outfile;
    outfile.open(filename);
    intT name;

    for (auto parent : gmap) {
        name = parent.first;

        for (auto child : parent.second) {
            outfile << name << " " << child << std::endl;
        }

    }
    outfile.close();
}

void convert_to_CSR(const GraphMap& gmap, CSR& cmat, int cmat_fmt) {

    intT n = gmap.size();
    cmat.fmt = cmat_fmt;
    cmat.rows = n;
    cmat.cols = n;
    cmat.nzcount = new intT[n];
    cmat.ja = new intT * [n];
    cmat.ma = new DataT * [n];

    //build the appropriate vectors from Mat;
    for (auto node : gmap) {

        intT parent = node.first; //parent node
        std::set<intT> tempset = node.second; //adiacency list for the node.
        intT nzs = tempset.size();
        cmat.nzcount[parent] = nzs;

        cmat.ja[parent] = new intT[nzs];
        cmat.ma[parent] = new DataT[nzs];

        std::copy(tempset.begin(), tempset.end(), cmat.ja[parent]);

        std::vector<DataT> temp_vec(tempset.size(), 1.);
        std::copy(temp_vec.begin(), temp_vec.end(), cmat.ma[parent]); //entries = 1s. GraphMap are unweighted (for now).
    }

}

int read_edgelist(std::string filename, CSR& cmat, int cmat_fmt, std::string delimiter = " ")
{
    std::ifstream infile;

    infile.open(filename);
    intT last_node = -1;
    intT current_node;
    std::string temp;
    std::vector<std::vector<intT>> holder = new std::vector < std::vector<intT>>();
    std::vector<intT> current_row;
    intT max_column = 0;
    
    while (infile.peek() == '#' or infile.peek() == '%') infile.ignore(2048, '\n');
    while (getline(infile, temp)) {

        int del_pos = temp.find(delimiter);
        int del_size = delimiter.length();

        std::string first_node_string = temp.substr(0, del_pos); //retrieve the part of the string before the delimiter
        current_node = stoi(first_node_string);
        temp.erase(0, del_pos + del_size);
        
        del_pos = temp.find(delimiter);
        std::string second_node_string = temp.substr(0, del_pos); //retrieve the part of the string after the delimiter
        child = stoi(second_node_string);
        max_column = std::max(max_column, child);

        if (current_node != last_node) current_row = new std::vector<intT>();
        if (current_node < last_node)
        {
            std::cout << "BAD FILE. CANNOT READ MATRIX" << std::endl;
        }
        last_node = current_node;

        current_row.push_back(child);
        holder.push_back(current_row);
    }

    cmat.fmt = cmat_fmt;
    intT main_dim = holder.size()
    intT second_dim = max_column + 1;
    cmat.rows = cmat_fmt ? second_dim : main_dim;
    cmat.cols = cmat_fmt ? main_dim : second_dim;
    cmat.nzcount = new intT[main_dim];
    cmat.ja = new intT * [main_dim];
    cmat.ma = new DataT * [main_dim];

    for (intT i = 0; i < holder.size(); i++)
    {
        auto row = holder[i];
        cmat.nzcount[i] = row.size();
        std::copy(row.begin(), row.end(), cmat.ja[i]);
        std::vector<DataT> temp_vec(row.size(), 1.);
        std::copy(temp_vec.begin(), temp_vec.end(), cmat.ma[i]); //entries = 1s. unweigthed
    }

    holder.clear();
    return 0;
}

void read_edgelist_DEPRECATED(std::string filename, CSR& cmat, int cmat_fmt, std::string delimiter = " ")
{
    GraphMap snap_graph;
    read_snap_format(snap_graph, filename, delimiter);         //Read into a GraphMap matrix from a .txt edgelist (snap format)
    MakeProper(snap_graph);
    convert_to_CSR(snap_graph, cmat, cmat_fmt);
}