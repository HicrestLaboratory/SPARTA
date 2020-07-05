#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <string>
#include <unistd.h>
#include <math.h>
#include <typeinfo>

// Utilities and system includes
#include <assert.h>
#include <helper_string.h>  // helper for shared functions common to CUDA Samples

// CUDA runtime
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>


// CUDA and CUBLAS functions
#include <helper_functions.h>
#include <helper_cuda.h>

#include "cuda_utilities.h"
#include "sparse_utilities.h"
#include "comp_mats.h"

using namespace std;

int main(int argc, char* argv[]) {


    if (typeid(DataT) != typeid(float)) {
        cout << "WARNING: only float supported for CUDA. Change DataT to float in sparse_utilities.h" << endl;
        return 1;
    }
    opterr = 0;

    int verbose = 3;

    int input_type = 4;
    int A_rows = 12;             //rows in the square input matrix;
    int A_cols = 8;
    int mat_A_fmt = 1;        //cuda needs column-major matrices
    int block_size = 4;     //block size for variable block matrix. Rows and columns must be evenly divisible by this;
    float sparsity = 0.5;   //sparsity of the input matrix;
    float block_sparsity = 0.5; //sparsity inside the blocks;
    string input_source;

    int B_cols = 5;    //number of columns in the output matrix;
    float B_sparsity = 0.8; //sparsity of the multiplication matrix

    float eps = 0.5;        //this value sets how different two rows in the same block can be.
                            //eps = 1 means only rows with equal structure are merged into a block
                            //eps = 0 means all rows are merged into a single block
    int seed = 123;




    srand(seed);

    int* A_row_part;
    int* A_col_part;




    //terminal options loop
    opterr = 0;
    char c;
    while ((c = getopt(argc, argv, "i:s:k:o:n:e:p:b:v:")) != -1)
        switch (c)
        {
        case 'i':// select input example
            input_type = stoi(optarg);
            //  1: Random CSR
            //  2: SNAP Edgelist
            //  3: MTX Format
            //  4: Random Variable Block matrix
            if ((input_type != 1) and (input_type != 4)) {
                input_type = 1;
                cout << "WARNING: CURRENTLY SUPPORTS ONLY i = 1 and i = 4. Using 1 (Random CSR)" << endl;
            }
            break;

        case 's': //select source file
            //has only effect for example 2 and 3;
            input_source = optarg;
            break;

        case 'k': //input matrix sparsity
            //has only effect for example 1 and 4
            sparsity = stof(optarg);
            if (sparsity < 0 or sparsity > 1) {
                fprintf(stderr, "Option -k tried to set sparsity outside of [0,1]");
                return 1;
            }
            break;

        case 'b': //block sparsity
        //has only effect for example 1 and 4
            block_sparsity = stof(optarg);
            if (block_sparsity < 0 or block_sparsity > 1) {
                fprintf(stderr, "Option -b tried to set block sparsity outside of [0,1]");
                return 1;
            }
            break;

        case 'n': //input matrix rows
             //has only effect for example 1 and 4
            A_rows = stoi(optarg);
            break;

        case 'o': //number of column of output matrix
            B_cols = stoi(optarg);
            break;

        case 'p': //size of blocks
            //ony used if i = 4, random VBS
            block_size = stoi(optarg);
            break;

        case 'v': //verbose
            verbose = stoi(optarg);
            break;

        case 'e': //epsilon used for matrix reordering;
            eps = stof(optarg);
            if (eps < 0. or eps > 1.) {
                fprintf(stderr, "Option -e tried to set epsilon outside of [0,1]");
                return 1;
            }
            break;

        case '?':
            fprintf(stderr, "Option -%c does not exists, or requires an argument.\n", optopt);
            return 1;
        default:
            abort();
        }


    //INPUT CONVERSION TO Compressed Sparse Row (CSR)

    CSR cmat_A; //this will hold the CSR matrix
    int cmat_A_fmt = 0;


    //INPUT EXAMPLE 1: RANDOM CSR
    //create a random sparse matrix
    if (input_type == 1) {
        DataT rand_mat[A_cols * A_rows];
        random_mat(rand_mat, A_rows, A_cols, sparsity); //generate random mat

        convert_to_CSR(rand_mat, A_rows, A_cols, mat_A_fmt, cmat_A, cmat_A_fmt);
        if (verbose > 0)
        {
            cout << "CREATED A RANDOM CSR" << endl;
        }
    }
    //______________________________________


    /*
    //NOT SUPPORTED
    //INPUT EXAMPLE 2: read graph in edgelist format into CSR
        if (input_type == 2){
            if (input_source.empty()) input_source = "testgraph.txt";

            read_snap_format(spmat, input_source);         //Read a CSR matrix from a .txt edgelist (snap format)
            cout << "IMPORTED A CSR FROM A SNAP EDGELIST" << endl;


        }
     //______________________________________
     */


     /*
     //NOT SUPPORTED
     //INPUT EXAMPLE 3: read from MTX format
         if (input_type == 3){
             //read from mtx
             if (input_source.empty()) input_source = "testmat.mtx";
             read_mtx_format(spmat, input_source); //read into CSR

             cout << "IMPORTED A CSR FROM MTX FILE" << endl;
             }


     //______________________________________
     */


     //INPUT EXAMPLE 4: create a random matrix with block structure
    if (input_type == 4) {

        //A_rows and sparsity have been previously set by options. Default: n = 20, sparsity = 0.5;

        DataT rand_block_mat[A_rows * A_cols];

        random_sparse_blocks_mat(rand_block_mat, A_rows, A_cols, mat_A_fmt, block_size, block_sparsity, sparsity);

        convert_to_CSR(rand_block_mat, A_rows, A_cols, mat_A_fmt, cmat_A, cmat_A_fmt);

        if (verbose > 0)
        {
            cout << "CREATED A RANDOM BLOCK MATRIX:"
                << " Rows = " << A_rows
                << " Columns: " << A_cols
                << " Block size: " << block_size
                << " Block sparsity: " << block_sparsity
                << " Entries sparsity: " << sparsity
                << endl;
        }
        //TODO optional: scramble the matrix row to see if the algo can reorder them.

    }

    //___________________________________________
    //*******************************************
    //		END OF INPUT
    //spmat must hold a proper CSR matrix at this point
    //******************************************

    //reorder the CSR matrix spmat and convert to a Block Sparse Matrix

        //TODO: reorder cmat through angle algorithm
        //      and retrieve row and block partition
        //      from the angle algorithm grouping;

    A_row_part = linspan(0, A_rows, block_size); //row and column partitions
    A_col_part = linspan(0, A_cols, block_size);


    //VBS matrix parameters
    int vbmat_blocks_fmt = 1;
    int vbmat_entries_fmt = 1; //cuda needs column-major matrices
    VBS vbmat_A;

    int block_rows = A_rows / block_size;
    int block_cols = A_cols / block_size;

    if ((A_rows % block_size != 0) or (A_cols % block_size != 0))
    {
        std::cout << "WARNING: The row or column dimension of the input matrix is not multiple of the block size " << std::endl;
    }

    convert_to_VBS(cmat_A,
        vbmat_A,
        block_rows, A_row_part,
        block_cols, A_col_part,
        vbmat_blocks_fmt, vbmat_entries_fmt);

    if (verbose > 0)
    {
        cout << "VBS matrix created:" << endl;
        matprint(vbmat_A);
    }

    //create a VBS with same structure as vbmat_A but which treats zero blocks as full blocks. Used for comparison.
    VBS vbmat_A_full;
    int no_zero_mode = 1;
    convert_to_VBS(cmat_A,
        vbmat_A_full,
        block_rows, A_row_part,
        block_cols, A_col_part,
        vbmat_blocks_fmt, vbmat_entries_fmt, no_zero_mode);

    if (verbose > 0)
    {
        cout << "VBS matrix (no zero blocks mode ON) created:" << endl;
        matprint(vbmat_A_full);
    }

    /*
    //*******************************************
    //        REPORT ON BLOCK STRUCTURE
    //******************************************
        ofstream CSV_out;
        CSV_out.open("output.txt");

        string CSV_header = "MatrixSize,OriginalSparsity,Divisions,NonzeroBlocks,AvgBlockHeight,AvgBHError,AvgBlockLength,AvgBLError,NonzeroAreaFraction,AverageBlockPopulation,ABPError,NewSparsity";
        CSV_out << CSV_header << endl;

        bool verbose = true; //print mat analysis on screen too?

        //TODO write this function
        //features_to_CSV(&vbmat, CSV_out, verbose);//write mat analysis on csv
        CSV_out.close();
    */


    //*******************************************
    //         MULTIPLICATION PHASE
    //___________________________________________
    //several ways of multiplying the sparse matrix
    //with a dense one, with benchmarks
    //******************************************

    //create a dense array matrix from spmat (for CUBLAS GEMM)
    DataT mat_A[A_rows * A_cols] = { 0 };

    convert_to_mat(cmat_A, mat_A, mat_A_fmt);
    if (verbose > 0)
    {
        std::cout << "Dense matrix A created:" << std::endl;
        matprint(mat_A, A_rows, A_cols, A_rows, mat_A_fmt);
    }

    //output format
    cout << fixed;

    if (verbose > 0)
    {
        cout << "\n \n **************************** \n STARTING THE MULTIPLICATION PHASE \n" << endl;
    }

    //creating a random matrix X
    int B_rows = A_cols;
    int mat_B_fmt = 1;

    DataT mat_B[B_rows * B_cols] = { 0 };
    random_mat(mat_B, B_rows, B_cols, B_sparsity);

    if (verbose > 0)        std::cout << "Random matrix B created:" << std::endl;
    if (verbose > 1)        matprint(mat_B, B_rows, B_cols, B_rows, mat_B_fmt);

    //creating the output matrix Y
	int C_rows = A_rows;
	int C_cols = B_cols;


    //--------------------------------------------
    //  dense-dense cublas gemm multiplication
    //--------------------------------------------

    DataT mat_Cgemm[C_rows * C_cols] = { 0 };
    int mat_Cgemm_fmt = 1;

    clock_t start_t = clock();

    cublas_gemm_custom (mat_A, A_rows, A_cols, A_rows, mat_B, B_cols, B_rows, mat_Cgemm, C_rows, 1.0f, 0.0f);

    double total_t = (clock() - start_t)/(double) CLOCKS_PER_SEC;
    
    if (verbose > 0)        cout << "Dense-Dense multiplication. Time taken: " << total_t << endl;
    if (verbose > 1)        
    {
        cout << "GEMM RESULT" << endl;  
        matprint(mat_Cgemm, C_rows, C_cols, C_rows, 1); 
    }

    //--------------------------------------------
    //      VBS x dense cublas multiplication	
    //--------------------------------------------

    DataT mat_Cblock[C_rows * C_cols];
    int mat_Cblock_fmt = 1;

    start_t = clock();

    cublas_blockmat_multiply(vbmat_A, mat_B, B_cols, B_rows, mat_Cblock, C_rows);

    total_t = (clock() - start_t) / (double)CLOCKS_PER_SEC;

    if (verbose > 0)
    {
        cout << "BlockSparse-Dense multiplication. Time taken: " << total_t << endl;
    }
    if (verbose > 1)
    {

        cout << "BLOCK RESULT" << endl;
        matprint(mat_Cblock, C_rows, C_cols, C_rows, 1);
    }

    int block_success = equal(C_rows, C_cols, mat_Cgemm, C_rows, mat_Cgemm_fmt, mat_Cblock, C_rows, mat_Cblock_fmt);
    if (block_success)
    {
        std::cout << "Block matrix multiplication test: SUCCESS" << std::endl;
    }
    else
    {
        std::cout << "Block matrix multiplication test: FAILED" << std::endl;
    }

    //--------------------------------------------
//      VBS x dense cublas multiplication	
//--------------------------------------------

    DataT mat_Cblock_full[C_rows * C_cols];
    int mat_Cblock_full_fmt = 1;

    start_t = clock();

    cublas_blockmat_multiply(vbmat_A_full, mat_B, B_cols, B_rows, mat_Cblock_full, C_rows);

    total_t = (clock() - start_t) / (double)CLOCKS_PER_SEC;

    if (verbose > 0)
    {
        cout << "BlockSparse-Dense multiplication (no zero mode ON). Time taken: " << total_t << endl;
    }
    if (verbose > 1)
    {

        cout << "BLOCK RESULT (no zero mode ON)" << endl;
        matprint(mat_Cblock_full, C_rows, C_cols, C_rows, 1);
    }

    int block_full_success = equal(C_rows, C_cols, mat_Cgemm, C_rows, mat_Cgemm_fmt, mat_Cblock_full, C_rows, mat_Cblock_full_fmt);
    if (block_full_success)
    {
        std::cout << "Block matrix multiplication test: SUCCESS" << std::endl;
    }
    else
    {
        std::cout << "Block matrix multiplication test: FAILED" << std::endl;
    }

}

//TODO CSR-dense cusparse multiplication

 
 