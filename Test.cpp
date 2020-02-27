#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <string>

#include "globheads.h"
#include "protos.h"
#include "utilities.h"
#include "mkl.h"


int main(int argc, char *argv[]) {
	
	double eps;
	SparMat spmat;
	if (argc > 1) {
		eps = stod(argv[1]);
	} 
	else {
		eps = 0.5; // 0 < eps < 1
	}
    //this value sets how different two rows in the same block can be.
    //eps = 1 means only rows with equal structure are merged into a block
    //eps = 0 means all rows are merged into a single block
    
	int rnd_dimension;
        if (argc > 2) {
                rnd_dimension = stod(argv[2]);
        }
        else {
                rnd_dimension = 10;
        }

    
    //create a random CSR matrix
    	Mat rand_mat;

   	float sparsity = 0.7;
	random_sparse_mat(rand_mat, rnd_dimension, sparsity); //generate random Mat
	convert_to_CSR(rand_mat, spmat);

    
//	cout << "READING GRAPH" << endl;    
    //Read a CSR matrix from a .txt edgelist (snap format)
//	read_snap_format(spmat, "testgraph.txt");
    
    
    //read from mtx
//	read_mtx_format(spmat, "testmat.mtx");


 //reorder the CSR matrix spamt and generate a Block Sparse Matrix
        VBSparMat vbmat;
        make_sparse_blocks(spmat, vbmat,eps);


//create a dense array matrix from spmat (for GEMM with MKL)
	Mat mat;
	int mat_n = spmat.n;
	double* mat_arr;
	mat_arr = new double[mat_n*mat_n];
	convert_from_CSR(spmat, mat);
	std::copy((mat.vec).begin(), (mat.vec).end(), mat_arr);

	cout << fixed;


//create a MKL sparse matrix from spmat
        sparse_matrix_t mkl_spmat;
        convert_to_MKL(spmat, mkl_spmat);

	

	
	cout << "PRINTING SPARSE MATRIX IN DENSE FORM" <<endl;
	matprint(mat_arr,mat_n,mat_n);    
	cout << "PRINTING SPARSE MATRIX IN BLOCK FORM" <<endl;    
	matprint(vbmat);

	ofstream CSV_out;
	CSV_out.open("output.txt");

	string CSV_header = "MatrixSize,OriginalSparsity,Divisions,NonzeroBlocks,AvgBlockHeight,AvgBHError,AvgBlockLength,AvgBLError,NonzeroAreaFraction,AverageBlockPopulation,ABPError,NewSparsity";
	CSV_out << CSV_header << endl;


	bool verbose = true; //print mat analysis
	features_to_CSV(&vbmat, CSV_out, verbose);//write mat analysis on csv
	CSV_out.close();
	

//MULTIPLICATION PHASE
	cout << "\n \n **************************** \n STARTING THE MULTIPLICATION PHASE \n" << endl; 
//creating the dense matrix X
	int X_rows = spmat.n;
	int X_cols = 1;	

	int k, seed = 4321;
  	srand(seed);
	double X[X_rows*X_cols];
  	for (k=0; k<X_rows*X_cols; k++) {
    		double x =  rand()%100;
    		X[k] = x/100;
  	}



//----------------------------
//creating the output matrix Y
	double Y_gemm[spmat.n * X_cols];
    	double Y_csr[spmat.n * X_cols];
    	double Y_block[spmat.n * X_cols] = {};
	double Y_batch[spmat.n * X_cols] = {};



//dense-dense mkl gemm multiplication
    
    clock_t start_t = clock();
    cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, mat_n, X_cols, mat_n, 1.0, mat_arr, mat_n, X, X_cols, 0, Y_gemm, X_cols);
    double total_t = (clock() - start_t)/(double) CLOCKS_PER_SEC;
        cout<<"Dense-Dense multiplication. Time taken: " << total_t<<endl;

    
//csr-dense mkl multiplication

	matrix_descr descr_spmat;
	descr_spmat.type = SPARSE_MATRIX_TYPE_GENERAL;
   	

	start_t = clock();
	mkl_sparse_d_mm (SPARSE_OPERATION_NON_TRANSPOSE, 1.0, mkl_spmat, descr_spmat, SPARSE_LAYOUT_ROW_MAJOR, X, X_cols, X_cols , 0.0, Y_csr, X_cols);
    	total_t = (clock() - start_t)/(double) CLOCKS_PER_SEC;
        cout<<"CSR-Dense multiplication. Time taken: " << total_t<<endl;


//column major version of X
	double X_c[X_rows*X_cols];

//	double X_test[X_cols*X_rows];

	convert_to_col_major(X, X_c, X_rows, X_cols);

//	convert_to_row_major(X_c,X_test,X_rows,X_cols);

//	matprint(X,X_rows,X_cols);
//	matprint(X_test,X_rows,X_cols);
//	cout << "test row/col conversions: " << are_equal(X,X_test, spmat.n*X_cols) << endl;

//------------------------------
	

//vbr-dense mkl multiplication	
	double Y_block_c[X_rows*X_cols] = {};

        start_t = clock();

	block_mat_multiply(vbmat, X_c, X_cols, Y_block_c);	
	
	total_t = (clock() - start_t)/(double) CLOCKS_PER_SEC;
	
	convert_to_row_major(Y_block_c,Y_block, spmat.n,X_cols);
	if(!are_equal(Y_block,Y_gemm, spmat.n*X_cols, 0.005)) cout << "VBR_DENSE Output check FAILED" << endl;
	
	cout <<"BlockSparse-Dense multiplication. Time taken: " << total_t<<endl;

//vbr-dense BATCH mkl multiplication
	double Y_batch_c[X_rows*X_cols] = {};

        start_t = clock();

        block_mat_batch_multiply(vbmat, X_c, X_cols, Y_block_c);

        total_t = (clock() - start_t)/(double) CLOCKS_PER_SEC;

        convert_to_row_major(Y_batch_c,Y_batch, spmat.n,X_cols);
        if(!are_equal(Y_block,Y_gemm, spmat.n*X_cols, 0.005)) cout << "VBR_BATCH Output check FAILED" << endl;

	cout <<"BlockSparse-Dense BATCH multiplication. Time taken: " << total_t<<endl;

//FINAL RESULTS
/*
	cout << "CSR RESULT" << endl;
        matprint(&Y_csr[0],spmat.n, X_cols);

	cout << "GEMM RESULT" << endl;
	matprint(&Y_gemm[0],spmat.n, X_cols);
	
	cout << "BLOCK RESULT" << endl;
	matprint(&Y_block[0],spmat.n, X_cols);

	cout << "BLOCK BATCH RESULT" << endl;
        matprint(&Y_batch[0],spmat.n, X_cols);
*/

}
