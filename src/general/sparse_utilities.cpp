#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <vector>
#include <algorithm>    // std::random_shuffle



#include "comp_mats.h"
#include "sparse_utilities.h"


typedef std::vector<int> svi;
typedef std::vector<DataT> svd;

// Matrix utilities

int IDX(int row, int col, int lead_dim, int fmt)
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

int is_empty(DataT* mat, int rows, int cols, int lead_dim, int fmt) {
    //check if a matrix is empty

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            if (mat[IDX(i, j, lead_dim, fmt)] != 0.) return 0;
        }
    }
    return 1;
}

int mat_cpy(DataT* in_mat, int in_rows, int in_cols, int in_lead_dim, int in_fmt, DataT* out_mat, int out_lead_dim, int out_fmt)
{
    //copy a matrix or submatrix in another matrix or submatrix of the same dimension. 

    int in_idx,out_idx;

    //TODO add error check

    for (int i = 0; i < in_rows; i++)
    {
        for (int j = 0; j < in_cols; j++)
        {
            in_idx  =   IDX(i, j, in_lead_dim, in_fmt);
            out_idx =   IDX(i, j, out_lead_dim, out_fmt);

            //TODO add out of bounds check
            out_mat[out_idx] = in_mat[in_idx];
        }
    }

    return 0;
}

int random_mat(DataT* mat, int rows, int cols, float sparsity)
{
    //create a random matrix with given sparsity and unitary entries;

    int nzs = (int)sparsity * rows * cols;
    svd entries = svd(rows * cols, 0);
    std::fill(entries.begin(), entries.begin() + nzs, 1);
    std::random_shuffle(entries.begin(), entries.end());
    std::copy(entries.begin(), entries.end(), mat);

}

int random_sparse_blocks_mat(DataT *mat, int rows, int cols, int fmt, int block_size, float block_sparsity, float block_entries_sparsity) 
{

    if ((rows % block_size != 0) or (cols % block_size != 0))
    {
        std::cout << "ERROR: matrix dimension must be a multiple of block size" << std::endl;
        return 1;
    }

    int n_blocks = rows * cols / (block_size * block_size);     //total number of blocks
    int nzblocks = (int)block_sparsity * n_blocks;              //total number of nonzero blocks

    svi blocks = svi(n_blocks, 0);              //will store 0 unless a block is nonzero;
    std::fill(blocks.begin(), blocks.begin() + nzblocks, 1);    //make nzblocks blocks nonzero;
    std::random_shuffle(blocks.begin(), blocks.end());          //put the nonzero blocks in random positions

    int mat_lead_dim = (fmt == 0) ? cols : rows;

    //put nonzerovalues in the Mat
    for (int ib = 0; ib < n_blocks; ib++) {//iterate through block rows
        for (int jb = 0; jb < n_blocks; jb++) { //iterate through block columns
            if (blocks[ib * n_blocks + jb] != 0) {
                //if block is nonempty, put random values in it;

                DataT tmp_block[block_size*block_size] = { 0 }; //temporary block
                random_mat(tmp_block, block_size, block_size, block_entries_sparsity); //make a random block with given sparsity

                int block_idx = IDX(ib * block_size, jb * block_size, mat_lead_dim, fmt); //find starting position of the block in mat
                mat_cpy(tmp_block, block_size, block_size, block_size, 0, mat + block_idx, mat_lead_dim, fmt); //copy entries from tmp_block to mat
            }
        }
    }


}

int matprint(DataT* mat, int rows, int cols, int lead_dim, int fmt)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            int idx = IDX(i, j, lead_dim, fmt);
            std::cout << mat[idx] << " ";
        }

        std::cout << std::endl;
    }
}

int sort_permutation(int* perm, int* arr, int n)
/*IN: perm, array of length n
      arr, array of length n
      n

  OUT: perm, the permutation that would sort arr
*/
{

    for (int i = 0; i < n; i++) {
        perm[i] = i; //sorted indices
    }
    //sort indices in perm by hash value
    std::sort(perm, perm + n,
        [&](const int& a, const int& b) {return (arr[a] < arr[b]); }
    );
}



//VBSfx utilities

int cleanVBS(VBSfx& vbmat)
{
    delete[] vbmat.nzcount;
    delete[] vbmat.jab;
    delete[] vbmat.mab;

    return 0;
}

int convert_to_VBSfx(DataT* mat, int mat_rows, int mat_cols, int mat_fmt, VBSfx& vbmat, int block_size, int vbmat_blocks_fmt, int vbmat_entries_fmt)
{

    //(For now only works with exact multiples.TODO add padding instead of removing elements)
    if ((mat_rows % block_size != 0) or (mat_cols % block_size != 0))
    {
        std::cout << "ERROR: matrix dimension must be a multiple of block size" << std::endl;
        return 1;
    }
    vbmat.rows = mat_rows/ block_size; //number of block-rows, round down.
    vbmat.cols = mat_cols/ block_size; //number of block-cols, round down.

    vbmat.blocks_fmt = vbmat_blocks_fmt;
    vbmat.entries_fmt = vbmat_entries_fmt;

    int vbmat_main_dim;
    int vbmat_compressed_dim;
    int i; //counter for main dimension
    int j; //counter for compressed dimension

    int* b_row_ptr; //pointer to current block_row
    int* b_col_ptr; //pointer to current block_column


    if (vbmat.blocks_fmt == 0)
    {
        //if Compressed Sparse Row
        vbmat_main_dim = vbmat.rows;
        vbmat_compressed_dim = vbmat.cols;

        // assign row and column pointers respectively
        // to counters for main and compressed dimension
        b_row_ptr = &i;
        b_col_ptr = &j;

    }
    else
    {
        //if Compressed Sparse Columns
        vbmat_main_dim = vbmat.cols;
        vbmat_compressed_dim = vbmat.rows;

        // assign column and row pointers respectively
        // to counters for main and compressed dimension
        b_row_ptr = &j;
        b_col_ptr = &i;
    }

    vbmat.nzcount = new int[vbmat_main_dim]; //initialize nzcount, which stores number of nonzero blocks for each block-row (-column)
    int mat_leading_dim = mat_fmt == 0 ? mat_cols : mat_rows;

    int matpos;
    svi jab;
    svd mab;

    //FIND BLOCK STRUCTURE--------------------------------------------------------------
    for (i = 0; i < vbmat_main_dim; i++ )        //loops trough main block dimension
    {
        for (j = 0; j < vbmat_compressed_dim; j++)     //loops through compressed block dimension
        {

            matpos = IDX((*b_row_ptr) * block_size, (*b_col_ptr) * block_size, mat_leading_dim, mat_fmt);    //find starting index of block in matrix
            if (!is_empty(mat + matpos, block_size, block_size, mat_leading_dim, mat_fmt))        //check if block is non-empty
            {
                vbmat.nzcount[i] += 1;  //one more nonzero block on the compressed dimension
                jab.push_back(j);       //store index of nonzero block
            }
        }
    }
    //----------------------------------------------------------------------------------

    vbmat.jab = new int[jab.size()];

    int total_nonzero_entries = jab.size() * block_size * block_size;
    vbmat.mab = new DataT[total_nonzero_entries];


    int mat_idx = 0; //keeps reading position for mat
    int vbmat_idx = 0; //keeps writing position for vbmat 
    int ja_count = 0; //keeps total nonzero blocks count;
    
    //COPY VALUES from mat to vbmat ------------------------------------------------------
    for (i = 0; i < vbmat_main_dim; i++)
    {
        for (int nzs = 0; nzs < vbmat.nzcount[i]; nzs++)
        {
            j = jab[ja_count];
            mat_idx = IDX((*b_row_ptr) * block_size, (*b_col_ptr) * block_size, mat_leading_dim, mat_fmt); //find starting index of block in matrix
            
            mat_cpy(mat + mat_idx, block_size, block_size, mat_leading_dim, mat_fmt, vbmat.mab + vbmat_idx, block_size, vbmat_entries_fmt); //write block from mat to vbmat.mab
            vbmat_idx += block_size * block_size;
            ja_count++;
        }
    }
    //------------------------------------------------------------------------------------

    return 0;
}

int convert_to_mat(const VBSfx& vbmat, DataT* out_mat, int out_mat_fmt)
{
    //input:
    //  vmat: a VBS matrix with fixed block dimension
    //  out_mat: an array of the proper dimension, filled with 0s; 
    int block_size = vbmat.block_size;

    //determine out_mat dimensions-------------------------
    int out_mat_rows = vbmat.rows * block_size;
    int out_mat_cols = vbmat.cols * block_size;
    int mat_leading_dim = out_mat_fmt == 0 ? out_mat_rows : out_mat_cols;
    //-----------------------------------------------------


    //determine vbmat dimensions-------------------------
    int vbmat_main_dim;
    int vbmat_compressed_dim;
    int i; //counter for main dimension
    int j; //counter for compressed dimension

    int* b_row_ptr; //pointer to current block_row
    int* b_col_ptr; //pointer to current block_column


    if (vbmat.blocks_fmt == 0)
    {
        //if Compressed Sparse Row
        vbmat_main_dim = vbmat.rows;
        vbmat_compressed_dim = vbmat.cols;

        // assign row and column pointers respectively
        // to counters for main and compressed dimension
        b_row_ptr = &i;
        b_col_ptr = &j;

    }
    else
    {
        //if Compressed Sparse Columns
        vbmat_main_dim = vbmat.cols;
        vbmat_compressed_dim = vbmat.rows;

        // assign column and row pointers respectively
        // to counters for main and compressed dimension
        b_row_ptr = &j;
        b_col_ptr = &i;
    }
    //--------------------------------------------------




    //FILL out_mat with the entries in vmat--------------
    int out_mat_idx = 0; //keeps read position for mat
    int vbmat_idx = 0; //keeps write position for vbmat 
    int ja_count = 0; //keeps total nonzero blocks count;

    for (i = 0; i < vbmat_main_dim; i++) //loop main dimension
    {
        for (int nzs = 0; nzs < vbmat.nzcount[i]; nzs++) //loop on compressed dimension
        {
            j = vbmat.jab[ja_count];
            out_mat_idx = IDX((*b_row_ptr) * block_size, (*b_col_ptr) * block_size, mat_leading_dim, out_mat_fmt); //find starting index of block in matrix

            mat_cpy(vbmat.mab + vbmat_idx, block_size, block_size, block_size, vbmat.entries_fmt, out_mat + out_mat_idx, mat_leading_dim, out_mat_fmt); //write block from vbmat.mab to mat
            vbmat_idx += block_size * block_size;
            ja_count++;
        }
    }
    //------------------------------------------------------

    return 0;
}

int convert_to_VBSfx(const CSR& cmat, VBSfx& vbmat, int block_size, int vbmat_block_fmt, int vbmat_entries_fmt)
{
    //WARNING: this does the additional step of converting to and from array. 
    //TODO: if necessary, make conversion efficient;

    int mat_rows = cmat.rows;
    int mat_cols = cmat.cols;
    int mat_size = mat_rows * mat_cols;
    int mat_fmt = 0;
    DataT mat[mat_size] = { 0 };
    convert_to_mat(cmat, mat, mat_fmt);
    convert_to_VBSfx(mat, mat_rows, mat_cols, mat_fmt, vbmat, block_size, vbmat_block_fmt, vbmat_entries_fmt);

    return 0;
}



//VBS utilities

int cleanVBS(VBS& vbmat)
{
    delete[] vbmat.nzcount;
    delete[] vbmat.jab;
    delete[] vbmat.mab;
    delete[] vbmat.rowpart;
    delete[] vbmat.colpart;

    return 0;
}

int convert_to_VBS(DataT* mat, int mat_rows, int mat_cols, int mat_fmt, VBS& vbmat, int block_rows, int* rowpart, int block_cols, int *colpart, int vbmat_blocks_fmt, int vbmat_entries_fmt)
{

    vbmat.rows = block_rows;
    vbmat.cols = block_cols;

    vbmat.blocks_fmt = vbmat_blocks_fmt;
    vbmat.entries_fmt = vbmat_entries_fmt;

    int vbmat_main_dim = vbmat.blocks_fmt == 0 ? vbmat.rows : vbmat.cols;
    int vbmat_compressed_dim = (vbmat.blocks_fmt == 0) ? vbmat.cols : vbmat.rows;
    int main_pos; //counter for main dimension
    int second_pos; //counter for compressed dimension

    vbmat.nzcount = new int[vbmat_main_dim]; //initialize nzcount, which stores number of nonzero blocks for each block-row (-column)
    int mat_leading_dim = mat_fmt == 0 ? mat_cols : mat_rows;

    int matpos;
    svi jab;
    svd mab;

    int main_pos, main_block_dim, second_pos, second_block_dim;
    int row, col, row_block_dim, col_block_dim;
    int total_nonzero_entries = 0;


    //FIND BLOCK STRUCTURE--------------------------------------------------------------
    for (i = 0; i < vbmat_main_dim; i++)        //loops trough main block dimension
    {
        main_pos = b_main_ptr[i];
        main_block_dim = b_main_ptr[i + 1] - main_pos;

        for (j = 0; j < vbmat_compressed_dim; j++)     //loops through compressed block dimension
        {

            second_pos = b_second_ptr[j];
            second_block_dim = b_second_ptr[j + 1] - second_pos;

            if (vbmat.block_fmt == 0)
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
            if (!is_empty(mat + matpos, row_block_dim, col_block_dim, mat_leading_dim, mat_fmt))        //check if block is non-empty
            {
                vbmat.nzcount[i] += 1;  //one more nonzero block on the compressed dimension
                tota_nonzero_entries += main_block_dim * second_block_dim;
                jab.push_back(j);       //store index of nonzero block
            }
        }
    }
    //----------------------------------------------------------------------------------

    vbmat.jab = new int[jab.size()];
    vbmat.mab = new DataT[total_nonzero_entries];


    int mat_idx = 0; //keeps reading position for mat
    int vbmat_idx = 0; //keeps writing position for vbmat 
    int jab_count = 0;

    //COPY VALUES from mat to vbmat ------------------------------------------------------
    for (i = 0; i < vbmat_main_dim; i++)
    {
        main_pos = b_main_ptr[i];
        main_block_dim = b_main_ptr[i + 1] - main_pos;
        for (int nzs = 0; nzs < vbmat.nzcount[i]; nzs++)
        {
            j = jab[jab_count];
            second_pos = b_second_ptr[j];
            second_block_dim = b_second_ptr[j + 1] - second_pos;

            if (vbmat.block_fmt == 0)
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

            block_leading_dim = (vbmat.entries_fmt == 0) ? second_block_dim : main_block_dim;
            
            mat_cpy(mat + mat_idx, row_block_dim, col_block_dim, mat_leading_dim, mat_fmt, vbmat.mab + vbmat_idx, block_leading_dim, vbmat_entries_fmt); //write block from mat to vbmat.mab
            vbmat_idx += main_block_dim*second_block_dim;
            jab_count++;
        }
    }
    //------------------------------------------------------------------------------------

    return 0;
}

int convert_to_mat(const VBS& vbmat, DataT* out_mat, int out_mat_fmt)
{
    //input:
    //  vmat: a VBS matrix
    //  out_mat: an array of the proper dimension, filled with 0s; 

    //determine out_mat dimensions-------------------------
    int out_mat_rows = vbmat.rowpart[vbmat.block_rows];
    int out_mat_cols = vbmat.colpart[vbmat.block_cols];
    int mat_leading_dim = out_mat_fmt == 0 ? out_mat_rows : out_mat_cols;
    //-----------------------------------------------------

    int vbmat_main_dim = vbmat.blocks_fmt == 0 ? vbmat.rows : vbmat.cols;
    int vbmat_compressed_dim = (vbmat.blocks_fmt == 0) ? vbmat.cols : vbmat.rows;
    int main_pos; //counter for main dimension
    int second_pos; //counter for compressed dimension

    int main_pos, main_block_dim, second_pos, second_block_dim;
    int row, col, row_block_dim, col_block_dim;

    int mat_idx = 0; //keeps writing position for mat
    int vbmat_idx = 0; //keeps writing position for vbmat 
    int ja_count = 0; //keeps total nonzero blocks count;

    //COPY VALUES from mat to vbmat ------------------------------------------------------
    for (i = 0; i < vbmat_main_dim; i++)
    {
        main_pos = b_main_ptr[i];
        main_block_dim = b_main_ptr[i + 1] - main_pos;
        for (int nzs = 0; nzs < vbmat.nzcount[i]; nzs++)
        {
            j = jab[nzs];
            second_pos = b_second_ptr[j];
            second_block_dim = b_second_ptr[j + 1] - second_pos;

            if (vbmat.block_fmt == 0)
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

            block_leading_dim = (vbmat.entries_fmt == 0) ? second_block_dim : main_block_dim;

            mat_cpy(vbmat.mab + vbmat_idx, row_block_dim, col_block_dim, block_leading_dim, vbmat_entries_fmt, mat + mat_idx, mat_leading_dim, mat_fmt); //write block from vbmat.mab to mat
            vbmat_idx += main_block_dim * second_block_dim; //update vbmat.mab reading index
        }
    }
    //------------------------------------------------------------------------------------

    return 0;

}

int convert_to_VBS(const CSR& cmat, VBS& vbmat, int block_rows, int* rowpart, int block_cols, int* colpart, int vbmat_block_fmt, int vbmat_entries_fmt)
{
    //WARNING: this does the additional step of converting to and from uncompressed array. 
    //TODO: if necessary, make conversion efficient;

    int mat_rows = cmat.rows;
    int mat_cols = cmat.cols;
    int mat_size = mat_rows * mat_cols;
    int mat_fmt = 0;
    DataT mat[mat_size] = { 0 };
    convert_to_mat(cmat, mat, mat_fmt);
    convert_to_VBS(mat, mat_rows, mat_cols, mat_fmt, vbmat, block_rows, rowpart, block_cols, colpart, vbmat_block_fmt, vbmat_entries_fmt);

    return 0;
}

int matprint(const VBS& vbmat)
{
    int mat_rows = vbmat.rowpart[vbmat.block_rows];
    int mat_cols = vbmat.colpart[vbmat.block_cols];

    DataT tmp_mat[mat_rows * mat_cols] = { 0 };
    convert_to_mat(vbmat, tmp_mat, 0);
    matprint(tmp_mat, mat_rows, mat_cols, cmat.cols, 0);
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


    int main_dim;
    main_dim = cmat.fmt == 0 ? cmat.rows : cmat.cols;

    for (int i = 0; i < main_dim; i++) {
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

    int main_dim;
    int compressed_dim;
    int i; //counter for main dimension
    int j; //counter for compressed dimension

    int* row_ptr; //pointer to current block_row
    int* col_ptr; //pointer to current block_column

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

    int out_lead_dim = out_mat_fmt == 0 ? cmat.cols : cmat.rows;
    
    //FILL THE OUTPUT MATRIX
    for (i = 0; i < main_dim; i++)
    {
        for (int nzs = 0; nzs < cmat.nzcount[i]; nzs++) 
        {
            j = cmat.ja[i][nzs]; //find column (row) index of next nonzero element
            DataT elem = cmat.ma[i][nzs]; //value of that element;

            int mat_idx = IDX(*row_ptr, *col_ptr, out_lead_dim, out_mat_fmt);
            out_mat[mat_idx] = elem;
        }
    }

    return 0;

}

int convert_to_CSR(const DataT* in_mat, int mat_rows, int mat_cols, int mat_fmt, CSR& cmat, int cmat_fmt)
{

    cmat.rows = mat_rows;
    cmat.cols = mat_cols;
    cmat.fmt = cmat_fmt;

    int cmat_main_dim;
    int cmat_compressed_dim;
    int i; //counter for main dimension
    int j; //counter for compressed dimension

    int* cmat_row_ptr; //pointer to current block_row
    int* cmat_col_ptr; //pointer to current block_column


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

    cmat.nzcount = new int[cmat_main_dim];
    cmat.ja = new int* [cmat_main_dim];
    cmat.ma = new DataT* [cmat_main_dim];

    int tmp_ja[cmat_compressed_dim];
    DataT tmp_ma[cmat_compressed_dim];

    int mat_lead_dim = mat_fmt == 0 ? cmat.cols : cmat.rows;

    for (i = 0; i < cmat_main_dim; i++) //loop through main dimension
    {
        int nzs = 0; //non-zero entries in this row (column) 
        for (j = 0; j < cmat_compressed_dim; j++) //scan compressed dimension for nonzero entries
        {
            int mat_idx = IDX(*cmat_row_ptr, *cmat_col_ptr, mat_lead_dim, mat_fmt);
            if (in_mat[mat_idx] != 0) //if the entry is nonzero add its idx to ja and its value to ma
            {
                tmp_ja[nzs] = j;
                tmp_ma[nzs] = in_mat[mat_idx];
                nzs++;
            }
        }

        cmat.nzcount[i] = nzs; 
        cmat.ja[i] = new int[nzs];
        cmat.ma[i] = new DataT[nzs];

        std::copy(tmp_ja, tmp_ja + nzs, cmat.ja[i]); 
        std::copy(tmp_ma, tmp_ma + nzs, cmat.ma[i]);

    }

    return 0;
}

int convert_to_CSR(const VBSfx& vbmat, CSR& cmat, int csr_fmt)
{
    //TODO: if necessary, make conversion efficient;

    int mat_rows = vbmat.rows * vbmat.block_size;
    int mat_cols = vbmat.cols * vbmat.block_size;
    int mat_size = mat_rows * mat_cols;
    int mat_fmt = 0;

    DataT mat[mat_size] = { 0 };
    convert_to_mat(vbmat, mat, mat_fmt);
    convert_to_CSR(mat, mat_rows, mat_cols, mat_fmt, cmat, csr_fmt);

    return 0;
}

int convert_to_CSR(const VBS& vbmat, CSR& cmat, int csr_fmt)
{
    //TODO: if necessary, make conversion efficient;

    int mat_rows = vbmat.rowpart[vbmat.block_rows];
    int mat_cols = vbmat.colpart[vbmat.block_cols];
    int mat_size = mat_rows * mat_cols;
    int mat_fmt = 0;

    DataT mat[mat_size] = { 0 };
    convert_to_mat(vbmat, mat, mat_fmt);
    convert_to_CSR(mat, mat_rows, mat_cols, mat_fmt, cmat, csr_fmt);

    return 0;
}

int matprint(const CSR& cmat)
{
    DataT tmp_mat[cmat.rows * cmat.cols] = { 0 };
    convert_to_mat(cmat, tmp_mat, 0);
    matprint(tmp_mat, cmat.rows, cmat.cols, cmat.cols, 0);
}

//TODO copyCSR

int transpose(const CSR& in_cmat, CSR& out_cmat, int new_fmt)
{
    int main_dim = in_cmat.fmt == 0 ? in_cmat.rows : in_cmat.cols; //main dimension of in_mat; secondary of out_mat; 
    int second_dim = in_cmat.fmt == 0 ? in_cmat.cols : in_cmat.cols; //secondary dimension of in_mat; main for out_mat;

    if (new_fmt != in_cmat.fmt)
    {
        out_cmat.fmt = new_fmt; //just change format instead of transposing
        out_cmat.rows = in_cmat.rows;
        out_cmat.cols = in_cmat.cols;

        //COPY CSR
        //RETURN
        std::cout << "FORMAT CHANGE NOT IMPLEMENTED YET. USE SAME FORMAT" << std::endl;
    }
    else
    {
        out_cmat.fmt = in_cmat.fmt;
        out_cmat.rows = in_cmat.cols;
        out_cmat.cols = in_cmat.rows;
    }


    out_cmat.nzcount = new int[second_dim];
    out_cmat.ja = new int* [second_dim];
    out_cmat.ma = new DataT *[second_dim];
    

    //find number of nonzero elements in each secondary row (which will be main row for the transpose); 
    for (int i = 0; i < main_dim; i++)
    {
        for (int nzs = 0; nzs < in_cmat.nzcount[i]; nzs++)
        {
            int j = in_cmat.ja[i][nzs]; //find column (row) index of next nonzero element
            out_cmat.nzcount[j] ++;
        }
    }

    int counter[second_dim] = { 0 };
   
    //initialize arrays in out_cmat
    for (int j = 0; j < second_dim; j++)
    {
        out_cmat.ja[j] = new int[out_cmat.nzcount[j]];
        out_cmat.ma[j] = new DataT[out_cmat.nzcount[j]];
    }


    int c = 0;
    for (int i = 0; i < main_dim; i++)
    {
        for (int nzs = 0; nzs < in_cmat.nzcount[i]; nzs++)
        {
            int j = in_cmat.ja[i][nzs]; //find in_cmat main_dim index of next nonzero element
            DataT elem = in_cmat.ma[i][nzs]; //value of that element;

            c = counter[j]; //progressively fill out_cmat main_dim
            out_cmat.ja[j][c] = i;
            out_cmat.ma[j][c] = elem;
        }
    }

}

int permute_CSR(CSR& cmat, int* perm, int dim) {
    //permutes rows (dim == 0), cols (dim == 1) or both (dim==2) of a matrix in CSR form;

    int main_dim = (cmat.fmt == 0) ? cmat.rows : cmat.cols;
    int second_dim = (cmat.fmt == 0) ? cmat.cols : cmat.rows;


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

        if (err_check != 0) return 1;
    }

    if (permute_second)
    {

        std::cout << "checksecond" << std::endl;

        int* ja;
        for (int i = 0; i < main_dim; i++)
        {
            ja = cmat.ja[i];
            int ja_len = cmat.nzcount[i];
            int idx_perm[second_dim];

            //change column indices to new values
            for (int j = 0; i < ja_len; j++)
            {
                ja[j] = perm[ja[j]];
                //index permutation FIX
            }

            int tmp_perm[ja_len];
            sort_permutation(tmp_perm, ja, ja_len); //find correct reorder of column indices.
            permute(cmat.ja[i], tmp_perm, ja_len); //sort ja[i]
            permute(cmat.ma[i], tmp_perm, ja_len); //permute ma[i] with the same permutation;
        }
    }
}

int hash_permute(CSR& cmat, int block_size, int* perm, int* group, int mode) 
{
    //finds a group structure and a permutation for the main dimension of a CSR mat

    // IN:
    //  cmat:        a matrix in CSR form
    //  block_size:  the number of elements to consider together when determining sparsity structure
    //              e.g if block_size = 8, every 8 element of secondary dimension will be considered nonzero if any of that is so
    //  mode:        0: at most one element per block is considered
    //              1: all elements in a block contribute to the hash
    // OUT:
    //  perm:        an array of length equal to cmat main dimension; stores a permutation of such dimension
    //  group:       an array of length equal to cmat main dimension; for each main row, stores the row group it belongs to


    int main_dim = cmat.fmt == 0 ? cmat.rows : cmat.cols;
    int second_dim = cmat.fmt == 0 ? cmat.cols : cmat.rows;

    if (block_size > second_dim)
    {
        std::cout << "ERROR. block_size must be smaller than matrix dimension" << std::endl;
        return 1;
    }

    int hashes[main_dim]; //will store hash values. The hash of a row (col) is the sum of the indices (mod block_size) of its nonzero entries

    for (int i = 0; i < main_dim; i++)
    {
        hashes[i] = hash(cmat.ja[i], cmat.nzcount[i], block_size, mode); //calculate hash value for each row
    }

    sort_permutation(perm, hashes, main_dim);

    int *ja_0, *ja_1;
    int len_0, len_1;
    int tmp_group = 0;
    group[0] = 0;

    for (int idx = 1; idx < main_dim; idx++) //scan main dimension in perm order and assign rows with same pattern to same group;
    {
        int i_0 = perm[idx - 1]; //previous element in perm order
        int i_1 = perm[idx]; //this element

        if (hashes[i_0] == hashes[i_1]) //only elements with the same hash have may have the same pattern
        {
            ja_0 = cmat.ja[i_0]; //previous row
            len_0 = cmat.nzcount[i_0]; //its length
            
            ja_1 = cmat.ja[i_1]; //this row
            len_1 = cmat.nzcount[i_1]; //its length

            if (check_same_pattern(ja_0, len_0, ja_1, len_1, block_size, mode)) //if new row has same pattern, put in the same group
            {
                group[idx] = tmp_group;
            }
            else                                                              // otherwise create new group
            {

                tmp_group += 1;
                group[idx] = tmp_group;

            }

        }
        else //when hash changes, create new group
        {
            tmp_group += 1;
            group[idx] = tmp_group;
        }
    }

}

int hash(int* arr, int a_len, int block_size, int mode)
{
    //evaluate hash function for an array of indices
    /* IN:
            arr: the array of indices.
            a_len: length fo the array;
            block_size: elements in the same block give the same contribution to hash
            mode:  0: at most one element per block contribute to hash
                   1: all elements in a block contribute to hash
        OUT: 
            int hash : the hash is the sum of the indices of nonzero blocks (indices are counted from 1, to avoid ignoring the 0 idx);
            */
       

    int nzs = 0;
    int tmp_idx = -1;
    int hash = 0;
    while (nzs < a_len)
    {
        int j = arr[nzs] / block_size;
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

int hash(int* arr, int a_len, int* block_partition, int blocks, int mode)
{
    //evaluate hash function for an array of indices
    /* IN:
            arr: the array of indices.
            a_len: length fo the array;
            block_partition: start position of block i; elements in the same block give the same contribution to hash
            blocks:  number of blocks;
            mode:  0: at most one element per block contribute to hash
                   1: all elements in a block contribute to hash
        OUT:
            int hash : the hash is the sum of the indices of nonzero blocks (indices are counted from 1, to avoid ignoring the 0 idx);
            */


    int nzs = 0;
    int hash = 0;
    
    int prev_idx = -1;
    int block_idx = 0;
    
    while (nzs < a_len)
    {
        while (arr[nzs] >= block_partition[bloc_idx + 1])
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

int check_same_pattern(int* arr0, int len0, int* arr1, int len1, int block_size, int mode)
{
    //check if two arrays of indices have the same pattern
    //PATTERN IS DEFINED BLOCKS OF FIXED DIMENSION block_size


    int i = 0;
    int j = 0;
    int b_idx0 = 0;
    int b_idx1 = 0;

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

int check_same_pattern(int* arr0, int len0, int* arr1, int len1, int* block_partition, int blocks, int mode)
{
    //check if two arrays of indices have the same pattern
    //PATTERN IS DIVIDED IN "blocks" BLOCKS OF VARIABLE DIMENSION, GIVEN BY block_partition; 

    int i = 0;
    int j = 0;

    int block_idx_0 = 0; //block_idx for arr0
    int block_idx_1 = 0; //block_idx for arr1


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
        if (b_idx0 != b_idx1)
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

//TODO check_has_pattern (int* arr, int len, int* pattern, int *block_partition, int_blocks, int mode);

int main()
{
    int rows = 10;
    int cols = 5;
    int fmt = 0;
    DataT mat[rows * cols] = { 0 };
    mat[1] = 1.;
    mat[15] = 2.;
    mat[28] = 3.;

    matprint(mat, rows, cols, cols, fmt);
    std::cout << "printed mat" << std::endl;

    CSR cmat; 
    int cmat_fmt = 0;
    convert_to_CSR(mat, rows, cols, fmt, cmat, cmat_fmt);

    matprint(cmat);
    std::cout << "mat converted to CSR" << std::endl;

    int perm1[cols] = { 1,2,3,0,4 };
    int perm2[rows] = { 1,2,0,9,8,7,3,4,6,5 };
    
    permute_CSR(cmat, perm2, 0);
    std::cout << "mat rows permuted" << std::endl;
    matprint(cmat);
    
    permute_CSR(cmat, perm1, 1);
    std::cout << "mat cols permuted" << std::endl;
    matprint(cmat);



}