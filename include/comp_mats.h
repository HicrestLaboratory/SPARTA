#pragma once
typedef float DataT; //precision for matrix entries
typedef long int intT;
typedef unsigned long int long_intT;


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
};

struct VBSfx {

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