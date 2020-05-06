#ifndef __COMPRESSED_MATRICES__
#define __COMPRESSED_MATRICES__

typedef struct CSR {
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
    int rows;      /* number of rows                                        */
    int cols;      /* number of cols                                        */
    int* nzcount;  /* number of nonzero entry in each row (column)          */
    int** ja;      /* pointer-to-pointer to store column (row) indices      */
    DataT** ma;   /* pointer-to-pointer to store nonzero entries           */
    int fmt;

};

typedef struct VBS {

    /*--------------------------------------------------------------
        Variable Block Sparse (VBS) matrix.
        The matrix doesn't need to be square,
        row and column partitions may differ.

        A VBS matrix has a main dimension and a compressed dimension.
        partition of rows and columns are stored respectively in row_part and col_part;
        specifically, row_part[i] holds the number of rows in the partition element i; 

        for each block row (or column) in the main dimension, we store:
            nzcount: the number of nonzero blocks
            jab: the indices of nonzero blocks, block row by block row (or column by column, see blocks_fmt)
            mab: the values of nonzero blocks, stored in row_major order (or column_major, see entries_fmt)
                 vectorized blocks are stored consecutively in the mab array.



------------------------------------------------------------                        */

    int block_side;         /*side lenght of blocks                                 */
    int rows;  	            /* the block row dimension of the matrix    	        */
    int cols;	            /* the block column dimension of the matrix   	        */
    int* nzcount;	        /* number of nonzero blocks in each block-row (-column) */
    int* jab;              /* block-columns (-rows) indices of nonzero blocks      */
    int* row_part;              /*cumulative number of row up to row partition element i*/
    int* col_part;              /*cumulative number of row up to row partition element i*/

    BData* mab;             /* array containing all entries, block by block	        */
    int entries_fmt;         /* storage format inside blocks:
                                0: row-major
                                1: column_major                                     */
    int blocks_fmt;          /* storage format of blocks:
                                    0: row-major
                                    1: column_major                                     */
};

typedef struct VBSfx {

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

    int block_side;         /*side lenght of blocks                                 */
    int rows;  	            /* the block row dimension of the matrix    	        */
    int cols;	            /* the block column dimension of the matrix   	        */
    int* nzcount;	        /* number of nonzero blocks in each block-row (-column) */
    int* jab;              /* block-columns (-rows) indices of nonzero blocks      */
    BData* mab;             /* array containing all entries, block by block	        */
    int entries_fmt;         /* storage format inside blocks:
                                0: row-major
                                1: column_major                                     */
    int blocks_fmt;          /* storage format of blocks:
                                    0: row-major
                                    1: column_major                                     */
};