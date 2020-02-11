#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "globheads.h"
#include "protos.h"
#include "utilities.h"
#include "mkl.h"

int main() {

	SparMat spmat;

	double eps = 0.5; // 0 < eps < 1
    //this value sets how different two rows in the same block can be.
    //eps = 1 means only rows with equal structure are merged into a block
    //eps = 0 means all rows are merged into a single block
    
    
    /*
    //create a random CSR matrix
    Mat mat;
    const int row_dimension = 500;
    float sparsity = 0.7;
	random_sparse_mat(mat, row_dimension, sparsity); //generate random Mat
	convert_to_CSR(mat, spmat);
    */
    
    /*
    //Read a CSR matrix from a .txt edgelist (snap format)
	read_snap_format(spmat, "testgraph.txt");
    */
    
    //read from mtx
    read_mtx_format(spmat, "testmat.mtx");

    sparse_matrix_t A;
    convert_to_MKL(spmat, A){

    //reorder the CSR matrix spmt and generate a Block Sparse Matrix
    VBSparMat vbmat;
    make_sparse_blocks(spmat, vbmat,eps);
    
    
	ofstream CSV_out;
	CSV_out.open("output.txt");

	string CSV_header = "MatrixSize,OriginalSparsity,Divisions,NonzeroBlocks,AvgBlockHeight,AvgBHError,AvgBlockLength,AvgBLError,NonzeroAreaFraction,AverageBlockPopulation,ABPError,NewSparsity";
	CSV_out << CSV_header << endl;
    
    bool verbose = true; //print mat analysis
	features_to_CSV(&vbmat, CSV_out, verbose);//write mat analysis on csv
	CSV_out.close();
}
