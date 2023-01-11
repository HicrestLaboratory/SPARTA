#pragma once


#include <string>
#include <vector>
#include <fstream>
#include "input.h"


typedef long int intT;
typedef float DataT; //precision for input matrices entries
typedef float DataT_C; //precision for result matrix entries
//types for the two multiplied matrices and the result matrix. available choices are: 
//          INT_8, INT_32
//          FLOAT_32, FLOAT_32

struct CSR
{
    /*--------------------------------------------------------------
    | Compressed sparse row (CSR) matrix format,
    | generally not square
    |==========
    |       cmat  = a CSR struct
    |       rows  = # of columns of the matrix
    |       cols  = # of columns of the matrix
            pattern_only   = 1: sparsity pattern only (boolean matrix)
    |               1: data and pattern
    |--------------------------------------------------------------------   */
    intT rows;      /* number of rows                                        */
    intT cols;      /* number of cols                                        */
    intT* nzcount;  /* number of nonzero entry in each row (column)          */
    intT** ja;      /* pointer-to-pointer to store column (row) indices      */
    DataT** ma;   /* pointer-to-pointer to store nonzero entries           */

    bool pattern_only; // 1 if the matrix is patter only; 0 otherwise

    void clean();
    void read_from_edgelist(std::ifstream& infile, std::string delimiter = " ", bool pattern_only = true);
    void reorder(std::vector<intT> grouping);
    void reorder_by_degree(bool descending = true);
    void permute_rows(std::vector<intT> permutation);
    void scramble();

    std::vector<intT> get_VBR_nzcount(const std::vector<intT> &grouping, intT block_col_size = 1);
    std::vector<intT> get_VBR_nzcount(const std::vector<intT> &row_partition, const std::vector<intT> &row_permutation, intT block_col_size = 1);
    void print(intT verbose = 0);
    intT nztot()
    {
        intT nztot = 0;
        for (intT i = 0; i < rows; i++)
        {
            nztot += nzcount[i];
        }
        return nztot;
    }


    //constructor for edgelist data
    CSR(std::ifstream& infile, std::string delimiter = " ", bool pattern_only = true)
    {
        read_from_edgelist(infile, delimiter, pattern_only);
    }


    //constructor from command line object
    CSR(CLineReader &cli)
    {
        std::ifstream fin;
        fin.open(cli.filename_);
        read_from_edgelist(fin, cli.reader_delimiter_, cli.pattern_only_);
        switch(cli.reorder_)
        {
            case 1:
                reorder_by_degree(true);
                break;
            case -1:
                reorder_by_degree(false);
                break;
            case 2:
                scramble();
                break;
        }
    }
    
    
    //destructor cleans all arrays
    ~CSR()
    {
        clean();
    }

};

struct VBR
{
    intT rows;
    intT cols;
    intT block_rows;  	            /* the block row dimension of the matrix    	        */
    intT block_cols;	            /* the block column dimension of the matrix   	        */
    intT* nzcount;	        /* i-th elements is the number of nonzero blocks in block-row i*/
    intT* jab;              /* block-columns indices of nonzero blocks      */
    intT* row_part;         /*row_part[i] is the index of the first row of the i-th block. Last element is number of rows for convenience  */
    DataT* mab;             /* array containing all entries, block by block	        */
    intT block_col_size;                  
    intT nztot;              /* total number of nonzero elements in mab*/


    //destructor cleans all arrays
    ~VBR()
    {
        clean();
    }

    void clean();
    void print(int verbose = 0);

    int partition_check(const std::vector<intT> &candidate_part);
    void fill_from_CSR(const CSR& cmat,const std::vector<intT> &row_partition, intT block_size);
    void fill_from_CSR_inplace(const CSR& cmat,const std::vector<intT> &grouping, intT block_size);
    DataT* get_block_start(intT row_block_idx);
};