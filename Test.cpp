#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <string>

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


    //create a MKL sparse matrix from spmat
    	sparse_matrix_t mkl_spmat;
    	convert_to_MKL(spmat, mkl_spmat);

    //create a dense array matrix from spmat
	Mat mat;
	int mat_n = spmat.n;
	float* mat_arr;
	mat_arr = new float[mat_n*mat_n];
	convert_from_CSR(spmat, mat);
	std::copy((mat.vec).begin(), (mat.vec).end(), mat_arr);


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


//MULTIPLICATION PHASE

//creating the dense matrix X
	int X_rows = spmat.n;
	int X_cols = spmat.n;

	int k, seed = 4321;
  	srand(seed);
	float X[X_rows*X_cols];
  	for (k=0; k<X_rows*X_cols; k++) {
    		float x = rand()%100;
    		X[k] = x/100;
  	}
//----------------------------
//creating the output matrix Y
	float Y_gemm[spmat.n * X_rows];
    float Y_csr[spmat.n * X_rows];
    float Y_block[spmat.n * X_rows];




//dense-dense mkl gemm multiplication
    
    clock_t start_t = clock();
    cblas_sgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, mat_n, mat_n, mat_n, 1.0, mat_arr, mat_n, X, mat_n, 0, Y_gemm,  mat_n);
    double total_t = (clock() - start_t)/(double) CLOCKS_PER_SEC;
        cout<<"Dense-Dense multiplication. Time taken: " << total_t<<endl;

    
//csr-dense mkl multiplication
	start_t = clock();

	matrix_descr descr_spmat;
	descr_spmat.type = SPARSE_MATRIX_TYPE_GENERAL;
    	
	mkl_sparse_s_mm (SPARSE_OPERATION_NON_TRANSPOSE, 1.0, mkl_spmat, descr_spmat, SPARSE_LAYOUT_ROW_MAJOR, X, X_rows, spmat.n , 0.0, Y_csr, spmat.n);
    total_t = (clock() - start_t)/(double) CLOCKS_PER_SEC;
        cout<<"CSR-Dense multiplication. Time taken: " << total_t<<endl;

//vbr-dense explicit multiplication

	



}
