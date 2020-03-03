#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <string>
#include <unistd.h>
#include <math.h>


#include "globheads.h"
#include "protos.h"
#include "utilities.h"
#include "mkl.h"


//TODO use float, not double matrices;
int main(int argc, char *argv[]) {

    opterr = 0;
    int input_type = 1;
    
    int n = 20; //rows in the square input matrix;
    
    int out_columns = 5; //number of columns in the output matrix;
    
    float sparsity = 0.5; //sparsity of the input matrix;
    
    string input_source;
    
    float eps = 0.8;
    //this value sets how different two rows in the same block can be.
    //eps = 1 means only rows with equal structure are merged into a block
    //eps = 0 means all rows are merged into a single block
   

    opterr = 0;
    char c;
    while ((c = getopt (argc, argv, "i:s:k:o:n:e:")) != -1)
      switch (c)
        {
        case 'i':// select input example
            input_type = stoi(optarg);
            //  1: Random CSR
            //  2: SNAP Edgelist
            //  3: MTX Format
            //  4: Random Variable Block matrix
            if(input_type < 1 or input_type > 4){
                input_type = 0;
                cout<<"WARNING: invalid input reference. Using 1 (Random CSR)"<<endl;
            }
            break;
        
        case 's': //select source file
            //has only effect for example 2 and 3;
            input_source = optarg;
            break;
        
        case 'k': //input matrix sparsity
            //has only effect for example 1 and 4
            sparsity = stof(optarg);
                if(sparsity < 0 or sparsity > 1){
                    fprintf (stderr, "Option -k tried to set sparsity outside of [0,1]");
                    return 1;
                }
          break;
                
        case 'n': //input matrix dimension
             //has only effect for example 1 and 4
            n = stoi(optarg);
	    break;
        
        case 'o': //number of column of output matrix
            out_columns = stoi(optarg);
            break;

        case 'e': //epsilon used for matrix reordering;
            eps = stof(optarg);
            if(eps < 0. or eps > 1.){
                fprintf (stderr, "Option -e tried to set epsilon outside of [0,1]");
                return 1;
            }
	    break;
                
        case '?':
            fprintf (stderr, "Option -%c does not exists, or requires an argument.\n", optopt);
            return 1;
        default:
          abort ();
	}
    

//INPUT SHOULD BE ALWAYS CONVERTED TO CSR BEFORE FURTHER MANIPULATION

	SparMat spmat; //this will hold the CSR matrix


//INPUT EXAMPLE 1: RANDOM CSR
//create a random sparse matrix
    if (input_type == 1){
        Mat rand_mat;
        random_sparse_mat(rand_mat, n, sparsity); //generate random Mat
        convert_to_CSR(rand_mat, spmat);
        cout << "CREATED A RANDOM CSR" << endl;

    }
//______________________________________

//TODO add import error notification
//INPUT EXAMPLE 2: read graph in edgelist format into CSR
    if (input_type == 2){
        if (input_source.empty()) input_source = "testgraph.txt";
        
        read_snap_format(spmat, input_source);         //Read a CSR matrix from a .txt edgelist (snap format)
        
        cout << "IMPORTED A CSR FROM A SNAP EDGELIST" << endl;


    }
 //______________________________________
        
        
//INPUT EXAMPLE 3: read from MTX format
    if (input_type == 3){
        //read from mtx
        if (input_source.empty()) input_source = "testmat.mtx";
        read_mtx_format(spmat, input_source);
        
        cout << "IMPORTED A CSR FROM MTX FILE" << endl;

        }


//______________________________________
//INPUT EXAMPLE 4: create a random matrix with block structure
    if (input_type == 4){

	int n_block = 3; //number of blocks
	float k_block = sqrt(sparsity); //percentage of non-zero blocks,must always be greater than sparsity

	Mat rnd_bmat;
	random_sparse_blocks_mat(rnd_bmat, n, n_block, k_block, sparsity);
        
	convert_to_CSR(rnd_bmat, spmat);
        
    cout << "CREATED A RANDOM BLOCK MATRIX" << endl;


	//optional: scramble the matrix?
    }
        
//___________________________________________
//*******************************************
//		END OF INPUT
//spmat must hold a SparMat matrix at this point
//******************************************

//reorder the CSR matrix spmat and generate a Block Sparse Matrix


        VBSparMat vbmat;
        make_sparse_blocks(spmat, vbmat,eps);	

	cout<<"CSR permuted. VBSparMat created"<<endl;

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

	

	
//	cout << "PRINTING SPARSE MATRIX IN DENSE FORM" <<endl;
//	matprint(mat_arr,mat_n,mat_n);    
//	cout << "PRINTING SPARSE MATRIX IN BLOCK FORM" <<endl;    
//	matprint(vbmat);



//*******************************************
//        REPORT ON BLOCK STRUCTURE
//*******************************************
	ofstream CSV_out;
	CSV_out.open("output.txt");

	string CSV_header = "MatrixSize,OriginalSparsity,Divisions,NonzeroBlocks,AvgBlockHeight,AvgBHError,AvgBlockLength,AvgBLError,NonzeroAreaFraction,AverageBlockPopulation,ABPError,NewSparsity";
	CSV_out << CSV_header << endl;


	bool verbose = true; //print mat analysis on screen?
	features_to_CSV(&vbmat, CSV_out, verbose);//write mat analysis on csv
	CSV_out.close();
	

//*******************************************
//         MULTIPLICATION PHASE
//___________________________________________
//various ways of multiplying the sparse matrix
//with a dense one, with benchmarks
//******************************************


//TODO put all matrices in column-major format
	cout << "\n \n **************************** \n STARTING THE MULTIPLICATION PHASE \n" << endl; 
//creating the dense matrix X
	int X_rows = spmat.n;
	int X_cols = out_columns;

	int seed = 123;
  	srand(seed);
	double X[X_rows*X_cols];
  	for (int k=0; k<X_rows*X_cols; k++) {
    		double x =  rand()%100;
    		X[k] = x/100;
  	}

        double X_c[X_rows*X_cols]; //column major version of X
        convert_to_col_major(X, X_c, X_rows, X_cols);


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
//------------------------------
	

//vbr-dense mkl multiplication	
	double Y_block_c[X_rows*X_cols] = {};

        start_t = clock();

	block_mat_multiply(vbmat, X_c, X_cols, Y_block_c);	
	
	total_t = (clock() - start_t)/(double) CLOCKS_PER_SEC;
	
	convert_to_row_major(Y_block_c,Y_block, spmat.n,X_cols);
	
	cout <<"BlockSparse-Dense multiplication. Time taken: " << total_t<<endl;


//BATCH MULTIPLUCATION NOT WORKING
    
//vbr-dense BATCH mkl multiplication
	double Y_batch_c[X_rows*X_cols] = {};

        start_t = clock();

        block_mat_batch_multiply(vbmat, X_c, X_cols, Y_block_c);

        total_t = (clock() - start_t)/(double) CLOCKS_PER_SEC;

        convert_to_row_major(Y_batch_c,Y_batch, spmat.n,X_cols);

	cout <<"BlockSparse-Dense BATCH multiplication. Time taken: " << total_t<<endl;
 
 
 
//PRINT RESULTING MATRICES
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
