#pragma once

#include <string>

typedef float DataT; //precision for input matrices entries
typedef float DataT_C; //precision for result matrix entries
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
    |--------------------------------------------------------------------   */
    intT rows;      /* number of rows                                        */
    intT cols;      /* number of cols                                        */
    intT* nzcount;  /* number of nonzero entry in each row (column)          */
    intT** ja;      /* pointer-to-pointer to store column (row) indices      */
    DataT** ma;   /* pointer-to-pointer to store nonzero entries           */
    int job;

    void clean();
    void read_from_edgelist(std::ifstream& infile, std::string delimiter, bool pattern_only);
    void reorder(std::vector<intT> permutation);
    void print();

};
