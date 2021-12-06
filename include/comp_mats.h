#pragma once
typedef float DataT; //precision for matrix entries
typedef float DataT_C; //precision for matrix entries
//types for the two multiplied matrices and the result matrix. available choices are: 
//          INT_8, INT_32
//          FLOAT_32, FLOAT_32

typedef long int intT;

struct CSR {
    /*--------------------------------------------------------------
    | Compressed sparse row (CSR) matrix format,
    | generally not square
    |==========
    |       cmat  = a CSR struct
    |       rows  = # of columns of the matrix
    |       cols  = # of columns of the matrix
            job   = 0: pattern only
    |               1: data and pattern
            fmt   = 0: compressed sparse rows
                    1: compressed sparse columns
    |--------------------------------------------------------------------   */
    intT rows;      /* number of rows                                        */
    intT cols;      /* number of cols                                        */
    intT* nzcount;  /* number of nonzero entry in each row (column)          */
    intT** ja;      /* pointer-to-pointer to store column (row) indices      */
    DataT** ma;   /* pointer-to-pointer to store nonzero entries           */
    int fmt;

};

struct VBS {

    /*--------------------------------------------------------------
        Variable Block Sparse (VBS) matrix.
        The matrix doesn't need to be square,
        row and column partitions may differ.

        A VBS matrix has a main dimension and a compressed dimension.
        partition of rows and columns are stored respectively in row_part and col_part;
        e.g., row_part[i] holds the row_index where partition element i starts;

        for each block row (or column) in the main dimension, we store:
            nzcount: the number of nonzero blocks
            jab: the indices of nonzero blocks, block row by block row (or column by column, see blocks_fmt)
            mab: the values of nonzero blocks, stored in row_major order (or column_major, see entries_fmt)
                 vectorized blocks are stored consecutively in the mab array.



------------------------------------------------------------                        */

    intT block_rows;  	            /* the block row dimension of the matrix    	        */
    intT block_cols;	            /* the block column dimension of the matrix   	        */
    intT* nzcount;	        /* number of nonzero blocks in each block-row (-column) */
    intT nztot;              /* total number of nonzero element in mab*/
    intT* jab;              /* block-columns (-rows) indices of nonzero blocks      */
    intT* col_part;              /*cumulative number of cols up to start of col partition element i (first element = 0, last element col_part[block_cols] is total number of cols) */
    intT* row_part;              /*cumulative number of rows up to start of row partition element i (first element = 0, last element row_part[block_rows] is total number of rows)*/

    DataT* mab;             /* array containing all entries, block by block	        */
    int entries_fmt;         /* storage format inside blocks:
                                0: row-major
                                1: column_major                                     */
    int blocks_fmt;          /* storage format of blocks:
                                    0: row-major
                                    1: column_major                                     */

    intT cols()
    {
        if (col_part) return col_part[block_cols];
        else return -1;
    }

    intT rows()
    {
        if (row_part) return row_part[block_rows];
        else return -1;
    }

    intT block_width(intT jb)
    {
        intT res = -1;
        if (jb < block_cols) res = (col_part[jb + 1] - col_part[jb]);
        return res;
    }

    intT block_height(intT ib)
    {
        intT res = -1;
        if (ib <= block_rows) res = (row_part[ib + 1] - row_part[ib]);
        return res;
    }

    intT main_dim()
    {
        return blocks_fmt ? block_cols : block_rows;
    }

    intT compressed_dim()
    {
        return blocks_fmt ? block_rows : block_cols;
    }

    intT* main_ptr()
    {
        return blocks_fmt ? col_part : row_part;
    }

    intT* second_ptr()
    {
        return blocks_fmt ? row_part : col_part;
    }

    float get_avg_height()
    {

        intT main_pos, main_block_dim, second_pos, second_block_dim;
        intT row, col, row_block_dim, col_block_dim;

        intT vbmat_main_dim = main_dim();
        intT* b_main_ptr = main_ptr();
        intT* b_second_ptr = second_ptr();

        float total_length = 0;
        intT nz_blocks = 0;
        //COPY VALUES from vbmat to mat ------------------------------------------------------
        for (intT i = 0; i < vbmat_main_dim; i++)
        {
            main_pos = b_main_ptr[i];
            main_block_dim = b_main_ptr[i + 1] - main_pos;

            for (intT nzs = 0; nzs < nzcount[i]; nzs++) //iterate for all nonzero block in row i
            {
                intT j = jab[nz_blocks]; //column of current non-zero block
                second_pos = b_second_ptr[j];
                second_block_dim = b_second_ptr[j + 1] - second_pos;
                nz_blocks++;
                if (blocks_fmt) total_length += second_block_dim;
                else total_length += main_block_dim;
            }
        }
        //------------------------------------------------------------------------------------

        return total_length/nz_blocks;
    }

    intT nz_blocks()
    {
        intT vbmat_main_dim = main_dim();
        intT nz_blocks = 0;
        for (intT i = 0; i < vbmat_main_dim; i++)
        {
            nz_blocks += nzcount[i];
        }
        return nz_blocks;
    }

};

struct ncVBS 
{
    /*
        Non-conform Variable Block Sparse matrix (ncVBS)

    */

    intT rows;                  //number of rows                     
    intT* col_part;             /*cumulative number of cols up to start of col partition element i (first element = 0, last element col_part[block_cols] is total number of cols) */
    intT* nzcount;	            /* len: block_cols; number of nonzero rows in each block-colums */
    intT block_cols;	        /* the block column dimension of the matrix   	        */
    intT** nzindex;              /* lenght: block_cols. Each vector contains the indices of the nzcount[i] nonzero rows in that column. */
    DataT** mab;                 /* lenght: block_cols. Each vector contains the nonzero elements in that column block, row by row. */
    DataT* mab_data;              /* length: sum of blocks dimensions*/

    intT cols()
    {
        if (col_part) return col_part[block_cols];
        else return -1;
    }

    intT tot_elems()
    {
        intT elems = 0;
        intT col_size;
        for (int jb = 0; jb < block_cols; jb++)
        {
            col_size = col_part[jb + 1] - col_part[jb];
            elems += col_size * nzcount[jb];
        }
        return elems;
    }

    intT block_width(intT jb)
    {
        return (col_part[jb + 1] - col_part[jb]);
    }

    intT elems_in_block(intT jb)
    {
        return  (col_part[jb + 1] - col_part[jb]) * nzcount[jb];
    }

    intT nz_elems_in_block(intT jb)
    {
        intT elems = block_width(jb) * nzcount[jb];
        intT nz = 0;
        for (int i = 0; i < elems; i++)
        {
            if (mab[jb][i] != 0)
            {
                nz++;
            }
        }
        return nz;
    }
};


struct VBSfx {


    //DEPRECATED
    /*--------------------------------------------------------------
        Block Sparse (VBS) matrix with square blocks of fixed dimension.
        The matrix itself doesn't need to be square, but must contain 
        an integer number of square blocks.

        A VBS matrix has a main dimension and a compressed dimension.
        for each block row (or column) in the main dimension, we store:
            nzcount: the number of nonzero blocks
            jab: the indices of nonzero blocks, block row by block row (or column by column, see blocks_fmt)
            mab: the values of nonzero blocks, stored in row_major order (or column_major, see entries_fmt)
                 vectorized blocks are stored consecutively in the mab array.



------------------------------------------------------------                        */

    intT block_size;         /*side lenght of blocks                                 */
    intT block_rows;  	            /* the block row dimension of the matrix    	        */
    intT block_cols;	            /* the block column dimension of the matrix   	        */
    intT* nzcount;	        /* number of nonzero blocks in each block-row (-column) */
    intT* jab;              /* block-columns (-rows) indices of nonzero blocks      */
    DataT* mab;             /* array containing all entries, block by block	        */
    int entries_fmt;         /* storage format inside blocks:
                                0: row-major
                                1: column_major                                     */
    int blocks_fmt;          /* storage format of blocks:
                                    0: row-major
                                    1: column_major                                     */
};