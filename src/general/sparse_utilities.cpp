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

// Matrix utilities

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

///TODO: make version with arbitrary secondary dimension partition, instead of fixed block side
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


//NOT TESTED YET
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
    partition(row_part, 0, block_rows, row_block_size); //row and block partition for creating the VBS
    partition(col_part, 0, block_cols, col_block_size);


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
    vbmat.mab = new DataT[nz_tot];
    vbmat.jab = new intT[nz_blocks];
    
    
    intT nz_in_block = std::ceil(entries_density * size_of_block);
    //saves the indices of nonzero blocks into jab; saves the number of nz blocks per row (or col) into vbmat.nzcount;
    intT b = 0;
    DataT* mab_idx = vbmat.mab; //the mab array will is filled up to this pointer
    for (intT i = 0; i < main_dim; i++)
    {
        intT nzcount = 0;
        for (intT j = 0; j < compressed_dim; j++)
        {
            if (blocks[b] != 0)
            {
                std::fill(mab_idx, mab_idx + nz_in_block, 1); //fill part of the block with ones
                std::random_shuffle(mab_idx, mab_idx + size_of_block); //shuffle the block
                vbmat.jab[b] == i; //save the index of the block in jab
                nzcount += 1; //keep the count of nonzero blocks
                mab_idx += size_of_block; // jump to next block on mab 
            }
            b++;
        }
        vbmat.nzcount[i] = nzcount;
    }

}

int matprint(DataT* mat, intT rows, intT cols, intT lead_dim, int fmt)
{
    for (intT i = 0; i < rows; i++)
    {
        for (intT j = 0; j < cols; j++)
        {
            intT idx = IDX(i, j, lead_dim, fmt);
            std::cout << mat[idx] << " ";
        }

        std::cout << std::endl;
    }
}

int matprint(DataT* mat, intT rows, intT* row_part, intT row_blocks, intT cols, intT* col_part, intT col_blocks, intT lead_dim,  int fmt)
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
                    std::cout << mat[idx] << " ";
                }

                std::cout << " ";

            }

            std::cout << std::endl;
        }

        std::cout << std::endl;
        
    }

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

int partition(intT* arr, intT start, intT end, intT step)
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

intT count_groups(intT* grp, intT grp_len)
{
    intT* perm = new intT[grp_len];
    sort_permutation(perm, grp, grp_len);
    intT last_grp = -1;
    intT groups = 0;

    for (intT idx = 0; idx < grp_len; idx++)
    {
        intT i = perm[idx];
        if (grp[i] != last_grp)
        {
            groups += 1;
            last_grp = grp[i];
        }
    }

    delete[] perm;
    return groups;

}

int grp_to_partition(intT* grp, intT grp_len, intT* partition)
{
    // IN: 
    //    grp: an array
    //    grp_len: lenght of grp
    // OUT: 
    //    partition: partition similar entries of grp together (sorted by entry)
    intT* perm = new intT[grp_len];
    sort_permutation(perm, grp, grp_len);
    intT last_grp = -1;
    intT group = 0;

    for (intT idx = 0; idx < grp_len; idx++)
    {
        intT i = perm[idx];
        if (grp[i] != last_grp)
        {
            partition[group] = idx;
            group += 1;
            last_grp = grp[i];
        }
        
    }
    partition[group] = grp_len;

    delete[] perm;
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
}

int convert_to_VBS(DataT* mat, intT mat_rows, intT mat_cols, int mat_fmt, VBS& vbmat, intT block_rows, intT* row_part, intT block_cols, intT *col_part, int vbmat_blocks_fmt, int vbmat_entries_fmt, int no_zero_mode)
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

    intT vbmat_main_dim, vbmat_compressed_dim;
    intT* b_main_ptr, * b_second_ptr;

    if (vbmat.blocks_fmt == 0)
    {
        //if Compressed Sparse Row
        vbmat_main_dim = vbmat.block_rows;
        vbmat_compressed_dim = vbmat.block_cols;

        // assign row and column pointers respectively
        // to counters for main and compressed dimension
        b_main_ptr = vbmat.row_part;
        b_second_ptr = vbmat.col_part;

    }
    else
    {
        //if Compressed Sparse Columns
        vbmat_main_dim = vbmat.block_cols;
        vbmat_compressed_dim = vbmat.block_rows;

        // assign column and row pointers respectively
        // to counters for main and compressed dimension
        b_main_ptr = vbmat.col_part;
        b_second_ptr = vbmat.row_part;
    }

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


//More efficient version: TODO TEST
int convert_to_VBS(const CSR& cmat, VBS& vbmat, intT block_rows, intT* row_part, intT block_cols, intT* col_part, int vbmat_block_fmt, int vbmat_entries_fmt)
{

    intT mat_rows = cmat.rows;
    intT mat_cols = cmat.cols;
    matprint(cmat);
    
    intT cmat_main_dim = cmat.fmt == 0 ? cmat.rows : cmat.cols;
    intT cmat_second_dim = cmat.fmt == 0 ? cmat.cols : cmat.rows;

    init_VBS(vbmat, block_rows, row_part, block_cols, col_part, vbmat_block_fmt, vbmat_entries_fmt);

    intT vbmat_main_dim = vbmat_block_fmt == 0 ? block_rows : block_cols;
    intT vbmat_second_dim = vbmat_block_fmt == 0 ? block_cols : block_rows;

    intT total_blocks = block_rows * block_cols;
    
    intT** blocks_bookmark = new intT*[vbmat_main_dim]; //will identify nonzero blocks;
    for (intT i = 0; i < vbmat_main_dim; i++)
    {
        blocks_bookmark[i] = new intT[vbmat_second_dim]{ -1 };
    }

    //pointers for proper iteration depending on fmt
    intT current_block_row = 0, current_block_col = 0;
    intT* current_main_pos = new intT();
    intT* current_second_pos = new intT();
    if (vbmat_block_fmt == 0)
    {
        intT* current_main_pos = &current_block_row;
        intT* current_second_pos = &current_block_col;
    }
    else
    {
        intT* current_main_pos = &current_block_col;
        intT* current_second_pos = &current_block_row;
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
            

            while (row > row_part[current_block_row])
            {
                current_block_row++;
            }

            while (col > col_part[current_block_col])
            {
                current_block_col++;
            }    

            //flag the bookmark position (nonzero block)
            blocks_bookmark[*current_main_pos][*current_second_pos] = -2;
            
        }

    }
    //DEBUG 
    int test = 0;
    std::cout << "test" << test++ << std::endl;


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
    std::cout << "test" << test++ << std::endl;

    //fill jab (nonzero blocks positions)
    total_nz_blocks = 0;
    for (intT ib = 0; ib < vbmat_main_dim; ib++)
    {
        vbmat.nzcount[ib] = 0;

        for (intT jb = 0; jb < vbmat_second_dim; jb++)
        {
            intT row = vbmat.blocks_fmt == 0 ? ib : jb;
            intT col = vbmat.blocks_fmt == 0 ? jb : ib;
            if (blocks_bookmark[ib][jb] != -1)
            {
                vbmat.jab[total_nz_blocks] = jb;
                total_nz_blocks += 1;
                vbmat.nzcount[ib] += 1;
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

            std::cout << "test i:" << i << " j: " << j << std::endl;

            while (row > row_part[current_block_row])
            {
                current_block_row++;
            }

            while (col > col_part[current_block_col])
            {
                current_block_col++;
            }

            //put the value in the correct block at the correct position

            intT block_start = blocks_bookmark[*current_main_pos][*current_second_pos];
            intT block_rows = row_part[current_block_row + 1] - row_part[current_block_row];
            intT block_cols = col_part[current_block_row + 1] - col_part[current_block_row];
            intT relative_row_pos = row - row_part[current_block_row];
            intT relative_col_pos = col - col_part[current_block_col];

            intT block_leading_dim = vbmat.entries_fmt == 0 ? block_rows : block_cols;
            intT block_idx = block_start + IDX(relative_row_pos, relative_col_pos, block_leading_dim, vbmat.entries_fmt);
            vbmat.mab[block_idx] = cmat.ja[i][nzs];
        }

    }
    std::cout << "test" << test++ << std::endl;

    for (intT i = 0; i < vbmat_main_dim; i++)
    {
        delete[] blocks_bookmark[i];
    }
    delete[] blocks_bookmark;
    return 0;
}


//NOT TESTED YET
//TO BE FINISHED
/*
int convert_to_VBS_direct(const CSR& cmat, VBS& vbmat, intT block_rows, intT* row_part, intT block_cols, intT* col_part, int vbmat_blocks_fmt, int vbmat_entries_fmt, int no_zero_mode)
{
    //convert from CSR to VBS

    intT nz_tot = count_nnz(cmat);
    
    cmat_main_dim = (cmat.fmt == 0)? cmat.rows : cmat.cols;
    vbmat_main_dim = (vbmat_blocks_fmt == 0) ? block_rows : block_cols;

    init_VBS(vbmat, block_rows, row_part, block_cols, col_part, vbmat_blocks_fmt, vbmat_entries_fmt);

    for (intT i = 0; i < cmat_main_dim; i++)
    {

    }
    init_VBS(vbmat, nz_tot, nz_blocks, block_rows, row_part, block_cols, col_part, vbmat_blocks_fmt, vbmat_entries_fmt);

}
*/

int matprint(const VBS& vbmat)
{

    intT rows = vbmat.row_part[vbmat.block_rows];    
    intT cols = vbmat.col_part[vbmat.block_cols];

    DataT* temp_mat = new DataT[rows * cols]{0};

    convert_to_mat(vbmat, temp_mat, 0);
    
    matprint(temp_mat, rows, vbmat.row_part, vbmat.block_rows, cols, vbmat.col_part, vbmat. block_cols, cols, 0);
    delete[] temp_mat;
}

//TODO: efficient conversion CSR <-> VBS


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

int convert_to_CSR(const VBS& vbmat, CSR& cmat, int csr_fmt)
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


void read_mtx_format(CSR& cmat, std::string infilename, int cmat_fmt) {
    std::ifstream file(infilename);
    intT rows, cols, num_lines;

    // Ignore comments headers
    while (file.peek() == '%') file.ignore(2048, '\n');

    // Read number of rows and columns
    file >> rows >> cols >> num_lines;

    std::vector<DataT> tmp_vec(rows * rows, 0.0);

    // fill the matrix with data
    for (intT l = 0; l < num_lines; l++)
    {
        DataT data;
        intT row, col;
        file >> row >> col >> data;
        tmp_vec[(row - 1) + (col - 1) * rows] = data;
    }

    DataT* tmp_mat = new DataT[rows * rows];
    std::copy(tmp_vec.begin(), tmp_vec.end(), tmp_mat);
    convert_to_CSR(tmp_mat, rows, cols, 0, cmat, cmat_fmt);
    delete[] tmp_mat;

}

//TODO: unify hash_permute, angle-permute and conversion to VBS

int hash_permute(CSR& cmat, intT* compressed_dim_partition, intT* perm, intT* group, int mode)
{
    //finds a group structure and a permutation for the main dimension of a CSR mat
    //NOTE: this does NOT automatically permute the CSR

    // IN:
    //  cmat:        a matrix in CSR form
    //  block_size:  the number of elements to consider together when determining sparsity structure
    //              e.g if block_size = 8, every 8 element of secondary dimension will be considered nonzero if any of that is so
    //  mode:        0: at most one element per block is considered
    //              1: all elements in a block contribute to the hash
    // OUT:
    //  perm:        an array of length equal to cmat main dimension; stores a permutation of such dimension that leaves groups together
    //  group:       an array of length equal to cmat main dimension; for each main row, stores the row group it belongs to

    intT main_dim = cmat.fmt == 0 ? cmat.rows : cmat.cols;
    intT second_dim = cmat.fmt == 0 ? cmat.cols : cmat.rows;

    intT* hashes = new intT[main_dim]; //will store hash values. The hash of a row (col) is the sum of the indices (mod block_size) of its nonzero entries

    for (intT i = 0; i < main_dim; i++)
    {
        group[i] = -1;

        hashes[i] = hash(cmat.ja[i], cmat.nzcount[i], compressed_dim_partition, mode); //calculate hash value for each row
    }

    sort_permutation(perm, hashes, main_dim); //find a permutation that sorts hashes


    intT *ja_0, *ja_1;
    intT len_0, len_1;
    intT tmp_group = -1;

    for (intT idx = 0; idx < main_dim; idx++) //scan main dimension in perm order and assign rows with same pattern to same group;
    {
        intT i = perm[idx]; //counter i refers to original order. Counter idx to permuted one. 

        if (group[i] == -1) //if row is still unassigned
        {

            tmp_group++; //create new group
            group[i] = tmp_group; //assign row to group
            
            ja_0 = cmat.ja[i]; //the row in compressed sparse format
            len_0 = cmat.nzcount[i];//the row length
            for (intT jdx = idx + 1; jdx < main_dim; jdx++)
            {
                intT j = perm[jdx]; //counter j refers to original order. Counter jdx to permuted one. 
                if (hashes[j] != hashes[i]) break; //only row with the same hashes must be compared, and they are consecutive in the perm order
                if (group[j] != -1) continue; //only unassigned row must be compared

                ja_1 = cmat.ja[j];
                len_1 = cmat.nzcount[j];

                if (check_same_pattern(ja_0, len_0, ja_1, len_1, compressed_dim_partition, mode))
                {
                    group[j] = tmp_group; //assign row j to the tmp_group
                }
            }
        }
    }

    delete[] hashes;
    sort_permutation(perm, group, main_dim); //stores in perm a permutation that sorts rows by group
}

intT hash(intT* arr, intT a_len, intT block_size, int mode)
{
    //evaluate hash function for a equally partitioned array of indices
    /* IN:
            arr: the array of indices.
            a_len: length fo the array;
            block_size: elements in the same block give the same contribution to hash
            mode:  0: at most one element per block contribute to hash
                   1: all elements in a block contribute to hash
        OUT: 
            intT hash : the hash is the sum of the indices of nonzero blocks (indices are counted from 1, to avoid ignoring the 0 idx);
            */
       

    intT nzs = 0;
    intT tmp_idx = -1;
    intT hash = 0;
    while (nzs < a_len)
    {
        intT j = arr[nzs] / block_size;
        nzs++;
        if ((j == tmp_idx) and (mode == 0)) //if mode is 0, only one j per block is considered in the hash sum;
        {
            continue;
        }

        hash += j + 1;
        tmp_idx = j;
    }
    return hash;
}

intT hash(intT* arr, intT a_len, intT* block_partition, int mode)
{
    //evaluate hash function for a arbitrarily partitioned array of indices
    /* IN:
            arr: the array of indices (a compressed row)
            a_len: length of the array;
            block_partition: start position of block i; elements in the same block give the same contribution to hash
            mode:  0: at most one element per block contribute to hash
                   1: all elements in a block contribute to hash
        OUT:
            intT hash : the hash is the sum of the indices of nonzero blocks (indices are counted from 1, to avoid ignoring the 0 idx);
            */


    intT nzs = 0;
    intT hash = 0;
    
    intT prev_idx = -1;
    intT block_idx = 0;
    
    while (nzs < a_len)
    {
        while (arr[nzs] >= block_partition[block_idx + 1])
        {
            block_idx++;
        };

        nzs++;
        if ((block_idx == prev_idx) and (mode == 0)) //if mode is 0, only one element per block is considered in the hash sum;
        {
            continue;
        }

        hash += block_idx + 1;
        prev_idx = block_idx;
    }
    return hash;
}

int check_same_pattern(intT* arr0, intT len0, intT* arr1, intT len1, intT block_size, int mode)
{
    //check if two arrays of indices have the same pattern
    //PATTERN IS DEFINED BLOCKS OF FIXED DIMENSION block_size


    intT i = 0;
    intT j = 0;
    intT b_idx0 = 0;
    intT b_idx1 = 0;

    while ((i < len0) and (j < len1))
    {
        b_idx0 = arr0[i] / block_size;
        b_idx1 = arr1[j] / block_size;
        if (b_idx0 != b_idx1)
        {
            return 0;
        }

        i++;
        j++; 

        if (mode == 0) //if mode=0, skip all entries in a block after the first one;
        {
            while ((i < len0) and (b_idx0 == arr0[i] / block_size))
            {
                i++;
            }
            while ((j < len1) and (b_idx1 == arr1[j] / block_size))
            {
                j++;
            }
        }

    }
    if ((i < len0) or (j < len1))
    {
        return 0;
    }
    else
    {
        return 1;
    }

}

int check_same_pattern(intT* arr0, intT len0, intT* arr1, intT len1, intT* block_partition, int mode)
{
    //check if two arrays of indices have the same pattern
    //PATTERN IS DIVIDED IN BLOCKS OF VARIABLE DIMENSION, GIVEN BY block_partition; 

    intT i = 0;
    intT j = 0;

    intT block_idx_0 = 0; //block_idx for arr0
    intT block_idx_1 = 0; //block_idx for arr1


    while ((i < len0) and (j < len1))
    {
        while (arr0[i] >= block_partition[block_idx_0 + 1])
        {
            block_idx_0++;
        }
        while (arr1[j] >= block_partition[block_idx_1 + 1])
        {
            block_idx_1++;
        }
        if (block_idx_0 != block_idx_1)
        {
            return 0;
        }

        i++;
        j++;

        if (mode == 0) //if mode=0, skip all entries in a block after the first one;
        {
            while ((i < len0) and (arr0[i] < block_partition[block_idx_0 + 1]))
            {
                i++;
            }
            while ((j < len1) and (arr1[j] < block_partition[block_idx_1 + 1]))
            {
                j++;
            }
        }

    }
    if ((i < len0) or (j < len1))
    {
        return 0;
    }
    else
    {
        return 1;
    }

}


//TODO get_pattern (and other pattern related functions) can be made more efficient by storing compressed blocks instead of an array of blocks
int get_pattern(intT* arr0, intT len0, intT* block_partition, intT* pattern, int mode)
{
    /*get the pattern of a compressed array w.r.t. a given block_partition
    IN:
        arr0: the compressed array, storing index of nonzero elements of the uncompressed array.
        len0: length of the compressed array, i.e. # of nonzero elements in the uncrompressed array;
        block_partition, bp: the partition. Block i has boundary [ bp[i], bp[i+1] ) ;
        mode:      0: at most one element per block contributes
                   1: all elements in a block contribute

    OUT:
        pattern: the pattern (has length of block partition - 1)
    */

    intT i = 0;
    intT block_idx = 0; //block_idx for pattern
    intT in_block = 0;

    while (i < len0)
    {
        while (arr0[i] >= block_partition[block_idx + 1]) //check if idx out of block; in that case, procede to next block.
        {
            pattern[block_idx] = in_block;
            in_block = 0;
            block_idx++;
        }

        in_block += 1;
        if ((mode == 0) and (in_block > 1))
        {
            in_block = 1;
        }

        pattern[block_idx] = in_block;
        i++;

    }
    return 0;
}

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


int angle_hash_method(CSR& cmat, float eps, intT* compressed_dim_partition, intT nB, VBS& vbmat, int vbmat_blocks_fmt, int vbmat_entries_fmt, int mode)
{
    //create a VBS reordering the main (uncompressed) dimension of a CSR matrix according to the angle+hash algorithm
    //do not change the original array

    intT rows = cmat.rows;
    intT cols = cmat.cols;
    intT main_dim = (cmat.fmt == 0) ? rows : cols;


    intT* hash_perm = new intT[main_dim];
    intT* hash_grp = new intT[main_dim];

    hash_permute(cmat, compressed_dim_partition, hash_perm, hash_grp, mode);
    
    intT* angle_perm = new intT[main_dim];
    intT* angle_grp = new intT[main_dim];

    angle_method(cmat, eps, compressed_dim_partition, nB, hash_perm, hash_grp, angle_grp, mode);

    sort_permutation(angle_perm, angle_grp, main_dim); //find a permutation that sorts groups

    intT angle_main_grps;
    angle_main_grps = count_groups(angle_grp, main_dim);

    intT* angle_main_part = new intT[angle_main_grps + 1];

    grp_to_partition(angle_grp, main_dim, angle_main_part);
    

    CSR cmat_cpy;
    copy(cmat, cmat_cpy);


    permute_CSR(cmat_cpy, angle_perm, cmat_cpy.fmt); //permute the tmp CSR

    intT* row_part;
    intT row_blocks;
    intT* col_part;
    intT col_blocks;

    if (cmat.fmt == 0)
    {
        row_part = angle_main_part;
        row_blocks = angle_main_grps;
        col_part = compressed_dim_partition;
        col_blocks = nB;
    }
    else
    {
        col_part = angle_main_part;
        col_blocks = angle_main_grps;
        row_part = compressed_dim_partition;
        row_blocks = nB;
    }

    convert_to_VBS(cmat_cpy,
        vbmat,
        angle_main_grps, angle_main_part,
        col_blocks, col_part,
        vbmat_blocks_fmt, vbmat_entries_fmt);

    //cleaning
    cleanCSR(cmat_cpy);
    delete[] hash_perm;
    delete[] hash_grp;
    delete[] angle_perm;
    delete[] angle_grp;
}


//TODO fix
/*
int symmetric_angle_hash_method(CSR& cmat, float eps, VBS& vbmat, intT vbmat_blocks_fmt, intT vbmat_entries_fmt)
{
    //create a VBS reordering the main (uncompressed) dimension of a CSR matrix according to the angle+hash algorithm
    //do not change the original array

    intT rows = cmat.rows;
    intT cols = cmat.cols;
    if (rows != cols)
    {
        std::cout << "WARNING! symmetric_angel_hash_method only accepts symmetric matrices" << std::endl;
        return 1;
    }

    intT n = rows;

    intT hash_perm[n];
    intT hash_grp[n];
    compressed_dim_partition = linspan(0, n, 1);

    hash_permute(cmat, compressed_dim_partition, hash_perm, hash_grp, 0);

    CSR cmat_cpy;
    copy(cmat, cmat_cpy);

    permute_CSR(cmat_cpy, hash_perm, 2); //permute both rows and columns of the matrix according to hash_perm;
    permute(hash_grp, hash_perm, n); //permute hash_grp to be aligned with the CSR
    
    intT hash_group_count = count_groups(hash_grp, n);
    intT hash_partition[hash_group_count];
    grp_to_partition(hash_grp, n, hash_partition); //extract row and block partition from grouping;

    hash_perm = linspace(0, n - 1, 1); //make angle_perm into identity permutation


    intT angle_perm[n];
    intT angle_grp[n];

    angle_method(cmat, eps, hash_partition, hash_group_count, hash_perm, hash_grp, angle_grp, mode);

    sort_permutation(angle_perm, angle_grp, n); //find a permutation that sorts groups


    intT angle_main_grps;
    angle_main_grps = count_groups(angle_grp, n);

    intT angle_main_part[angle_main_grps];

    grp_to_partition(angle_grp, n, angle_main_part);

    permute_CSR(cmat_cpy, angle_perm, 2); //permute the tmp CSR

    intT* row_part;
    intT row_blocks;
    intT* col_part;
    intT col_blocks;

    if (cmat.fmt == 0)
    {
        row_part = angle_main_part;
        row_blocks = angle_main_grps;
        col_part = compressed_dim_partition;
        col_blocks = nB;
    }
    else
    {
        col_part = angle_main_part;
        col_blocks = angle_main_grps;
        row_part = compressed_dim_partition;
        row_blocks = nB;
    }

    convert_to_VBS(cmat_cpy,
        vbmat,
        angle_main_grps, angle_main_part,
        col_blocks, col_part,
        vbmat_blocks_fmt, vbmat_entries_fmt);

    cleanCSR(cmat_cpy);
}
*/

int angle_method(CSR& cmat, float eps, intT* compressed_dim_partition, intT nB,intT* in_perm, intT* in_group, intT *out_group,  int mode)
{
    /*
    COMPUTES A GROUPING AND PERMUTATION FOR A CSR MATRIX (ALONG ITS MAIN DIMENSION) 
    GIVEN A STARTING PERMUTATION AND A GROUPING; USES ANGLE METHOD WITH A FIXED PARTITION OF THE COMPRESSED DIMENSION.

    IN: 
        cmat: the CSR matrix
        compressed_dim_partition: the partition of the matrix along its compressed dimension. Element i is the start position of block i; last element is the length of the compressed dimension;
        nB: the number of blocks in the partition
        in_perm: a permutation of the main dimension (rows in the same group should be adjacent in this permutation); 
        in_group: a grouping along the main dimension;
        eps: rows with cosine greater than eps will be merged.
        mode: 0: consider at most one element per block when evaluating sparsity patterns
              1: consider all the elements in a block (but not their order) when evaluating sparsity patterns

    OUT: out_group will store a new grouping (compatible with in_group)
    */

    intT main_dim = cmat.fmt == 0 ? cmat.rows : cmat.cols;
    intT second_dim = cmat.fmt == 0 ? cmat.cols : cmat.rows;

    for (intT i = 0; i < main_dim; i++)//initialize out_group
    {
        out_group[i] = -1;
    }

    intT idx, jdx;
    intT this_group = -1; //the (out_)group the current row is in. 
    intT in_this_grp;
    intT* this_pattern = new intT[nB]{ 0 }; //the pattern of the current row or group of rows

    intT that_group;//the (in_)group the compared row is in
    intT in_that_grp;
    intT* that_pattern = new intT[nB]{ 0 };

    intT i, j;

    for (intT idx = 0; idx < main_dim; idx++) //Loop through (groups of) rows. Each one is confronted with all the unpaired ones to find those that will be merged.
    {
        i = in_perm[idx];           //idx counts in the permuted order. i counts in the original order;

        if (out_group[i] == -1)     //only consider still ungrouped rows;
        {
            this_group++;               //create new group

            intT* arr0 = cmat.ja[i];
            intT len0 = cmat.nzcount[i];
           
            jdx = idx; //idx to compare the row with. will only compare with elements after current row;
 
            in_this_grp = 0;
            while ((jdx < main_dim) and (in_group[in_perm[jdx]] == in_group[i])) //check for elements in the same in_group of i (they are consecutive in in_perm)
            {
                out_group[in_perm[jdx]] = this_group;   //assign elements in the same in_group to the same out_group;
                jdx++;
                in_this_grp++; //keep count of the elements in the group
            }

            get_pattern(arr0, len0, compressed_dim_partition, this_pattern, mode); //get the row pattern (stores into this_pattern)




            if (mode == 1)
            {
                for (intT b = 0; b < nB; b++)
                {
                    this_pattern[b] *= in_this_grp; //multiply entries of the pattern with number of rows with the same pattern (same in_group)
                }
            }
            intT norm_0 = norm2(this_pattern, nB); //squared norm of the pattern

            while(jdx < main_dim) //loop through not-analyzed rows to be paired with the group
            {
                j = in_perm[jdx];
                that_group = in_group[j];
                bool merge = false;

                if (out_group[j] == -1)         //only consider ungrouped rows;
                {
                    intT* arr1 = cmat.ja[j];
                    intT len1 = cmat.nzcount[j];
                    that_group = in_group[j];

                    get_pattern(arr1, len1, compressed_dim_partition, that_pattern, mode); //get the row pattern (store into that_pattern)

                    intT norm_1 = norm2(that_pattern, nB); //get norm of the pattern

                    float scal = scalar_product(this_pattern, nB, that_pattern);
                    if (scal*scal > eps * norm_0 * norm_1) //if cosine is > than epsilon, allow merge of groups
                    {
                        merge = true;
                    }
                }

                in_that_grp = 0; 
                while ((jdx < main_dim) and (in_group[in_perm[jdx]] == that_group)) //iterate over elements in the same group of the analyzed row j (j included)
                {
                    in_that_grp++;
                    if (merge) out_group[in_perm[jdx]] = this_group; //if merge was decided, add group to current group 
                    jdx++;
                }

                if (merge)
                {
                    for (intT b = 0; b < nB; b++)
                    {
                        if (mode == 0)
                        {
                            this_pattern[b] = std::max(this_pattern[b], that_pattern[b]);
                        }
                        else
                        {
                            this_pattern[b] += in_that_grp * that_pattern[b]; //update pattern with entries of merged rows
                        }
                    }
                    norm_0 = norm2(this_pattern, nB); //update norm of this pattern
                }
            }
        }
    }

    delete[] this_pattern;
    delete[] that_pattern;


}


void read_snap_format(GraphMap& gmap, std::string filename)
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
    while ((infile.peek() == '%') or (infile.peek() == '#')) infile.ignore(2048, '\n');

    //TODO: check import success
    //read the source node (row) of a new line
    while (getline(infile, temp, ' ')) {
        current_node = stoi(temp);

        if (gmap.count(current_node) == 0) { //new source node encountered
            gmap[current_node] = emptyset;
        }

        //get target node (column)
        getline(infile, temp);
        child = stoi(temp);
        gmap[current_node].insert(child);
    }
}

void read_snap_format(GraphMap& gmap, std::string filename, std::string delimiter)
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
    //read the source node (row) of a new line
    while (getline(infile, temp)) {

        intT del_pos = temp.find(delimiter);
        intT del_size = delimiter.size();

        std::string first_node_string = temp.substr(0, del_pos); //retrieve the part of the string before the delimiter
        current_node = stoi(first_node_string);

        if (gmap.count(current_node) == 0) { //new source node encountered
            gmap[current_node] = emptyset;
        }

        //get target node (column)
        std::string second_node_string = temp.substr(del_pos + del_size); //retrieve the part of the string after the delimiter
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
    cmat.ja = new intT* [n];
    cmat.ma = new DataT * [n];

    //build the appropriate vectors from Mat;
    for (auto node : gmap) {

        intT parent = node.first; //parent node
        std::set<intT> tempset = node.second; //adiacency list for the node.
        intT nzs = tempset.size();
        cmat.nzcount[parent] = nzs;

        cmat.ja[parent] = new intT [nzs];
        cmat.ma[parent] = new DataT [nzs];

        std::copy(tempset.begin(), tempset.end(), cmat.ja[parent]);

        std::vector<DataT> temp_vec(tempset.size(), 1.);
        std::copy(temp_vec.begin(), temp_vec.end(), cmat.ma[parent]); //entries = 1s. GraphMap are unweighted (for now).
    }

}