#pragma once


#include <string>
#include <vector>
#include <fstream>

typedef long int intT;
typedef float DataT; //precision for input matrices entries
typedef float DataT_C; //precision for result matrix entries
//types for the two multiplied matrices and the result matrix. available choices are: 
//          INT_8, INT_32
//          FLOAT_32, FLOAT_32

struct CSR {
    /*--------------------------------------------------------------
    | Compressed sparse row (CSR) matrix format,
    | generally not square
    |==========
    |       cmat  = a CSR struct
    |       rows  = # of columns of the matrix
    |       cols  = # of columns of the matrix
            job   = 0: sparsity pattern only (boolean matrix)
    |               1: data and pattern
    |--------------------------------------------------------------------   */
    intT rows;      /* number of rows                                        */
    intT cols;      /* number of cols                                        */
    intT* nzcount;  /* number of nonzero entry in each row (column)          */
    intT** ja;      /* pointer-to-pointer to store column (row) indices      */
    DataT** ma;   /* pointer-to-pointer to store nonzero entries           */
    
    //arrays that contains all the values
    DataT* ma_full;      //hosts the values of ma
    intT* ja_full;  // hosts the values of ja

    int job; // 0 if the matrix is patter only; 1 otherwise

    void clean();
    void read_from_edgelist(std::ifstream& infile, std::string delimiter, bool pattern_only);
    void reorder(std::vector<intT> permutation);
    void print();

};
