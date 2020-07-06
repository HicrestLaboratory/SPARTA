#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <string>
#include <unistd.h>
#include <math.h>
#include <typeinfo>
#include <iterator>
#include <algorithm>

// Utilities and system include
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
typedef std::vector<float> vec_d;
typedef std::vector<string> vec_str;
typedef std::vector<int> vec_i;

float mean(vec_d v)
{
    float m = 0.;
    for (auto t : v)
    {
        m += t;
    }
    m /= v.size();
    return m;
}

float std_dev(vec_d v)
{
    float m = mean(v);
    float s = 0.;
    for (auto t : v)
    {
        s += (t - m) * (t - m);
    }
    s = sqrt(s / v.size());
}

template <class myType>
int output_couple(string& names, string& values, string name, myType value)
{
    names += name + " ";
    values += to_string(value) + " ";
}

int main(int argc, char* argv[]) {


    if (typeid(DataT) != typeid(float)) {
        cout << "WARNING: only float supported for CUDA. Change DataT to float in sparse_utilities.h" << endl;
        return 1;
    }
    opterr = 0;

    int verbose = 3;

    int input_type = 4;
    int algos = -1;
    int A_rows = 12;             //rows in the square input matrix;
    int A_cols = 8;
    int mat_A_fmt = 1;        //cuda needs column-major matrices
    int block_size = 4;     //block size for variable block matrix. Rows and columns must be evenly divisible by this;
    float sparsity = 0.5;   //sparsity of the input matrix;
    float block_sparsity = 0.5; //sparsity inside the blocks;
    string input_source;
    int scramble = 1; //scramble the input matrix?

    int B_cols = 5;    //number of columns in the output matrix;
    float B_sparsity = 1.; //sparsity of the multiplication matrix

    float eps = 0.5;        //this value sets how different two rows in the same block can be.
                            //eps = 1 means only rows with equal structure are merged into a block
                            //eps = 0 means all rows are merged into a single block
    int seed = 123;
    float precision = 0.0001;        //precision for float equality check

    int warmup = 0;         //number of warmup experiments
    int experiment_reps = 5; //number of non-warmup repetitions
    int algo = -1;           //algorithm choice (-1: all)
    int correct_check = 0;




    srand(seed);

    int* A_row_part;
    int* A_col_part;




    //terminal options loop
    opterr = 0;
    char c;
    while ((c = getopt(argc, argv, "i:s:q:e:p:m:n:k:b:v:")) != -1)
        switch (c)
        {
        case 'i':// select input example
            input_type = stoi(optarg);
            //  1: Random CSR
            //  2: SNAP Edgelist
            //  3: MTX Format
            //  4: Random Variable Block matrix
            if ((input_type != 1) and (input_type != 4) and (input_type != 3)) {
                input_type = 4;
                cout << "WARNING: CURRENTLY SUPPORTS ONLY i = 1,3,4. Using 4 (Random CSR)" << endl;
            }
            break;

        case 's': //select source file
            //has only effect for example 2 and 3;
            input_source = optarg;
            break;

        case 'q': //input matrix sparsity
            //has only effect for example 1 and 4
            sparsity = stof(optarg);
            if (sparsity < 0 or sparsity > 1) {
                fprintf(stderr, "Option -k tried to set sparsity outside of [0,1]");
                return 1;
            }
            break;

        case 'b': //block sparsity
        //has only effect for example 4
            block_sparsity = stof(optarg);
            if (block_sparsity < 0 or block_sparsity > 1) {
                fprintf(stderr, "Option -b tried to set block sparsity outside of [0,1]");
                return 1;
            }
            break;

        case 'm': //input matrix rows
             //has only effect for example 1 and 4
            A_rows = stoi(optarg);
            break;

        case 'n': //input matrix rows
     //has only effect for example 1 and 4
            B_cols = stoi(optarg);
            break;

        case 'k': //input matrix rows
 //has only effect for example 1 and 4
            A_cols = stoi(optarg);
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

    //TODO -h HELP

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


     //INPUT EXAMPLE 3: read from MTX format
    if (input_type == 3) 
    {
        //read from mtx
        if (input_source.empty()) input_source = "testmat.mtx";
        read_mtx_format(cmat_A, input_source, cmat_A_fmt); //read into CSR

        if (verbose > 0)
        {
            cout << "IMPORTED A CSR FROM MTX FILE" << endl;
        }
    }



     //______________________________________


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


    //scramble the original matrix
    if (scramble)
    {
        int* random_permutation = randperm(A_rows);
        permute_CSR(cmat_A, random_permutation, 0);
    }

    //*******************************************
    //	 CREATES (several) VBS from the CSR, EXTRACT FEATURES
    //******************************************

    string output_names;
    string output_values;


    float A_sparsity = ((float) count_nnz(cmat_A))/(A_rows*A_cols);


    output_couple(output_names, output_values, "input_type", input_type);
    output_couple(output_names, output_values, "A_rows", A_rows);
    output_couple(output_names, output_values, "A_cols", A_cols);
    output_couple(output_names, output_values, "A_entries_sparsity", A_sparsity);
    output_couple(output_names, output_values, "A_block_sparsity", block_sparsity);
    output_couple(output_names, output_values, "A_block_size", block_size);

    output_couple(output_names, output_values, "B_cols", B_cols);
    output_couple(output_names, output_values, "B_sparsity", B_sparsity);


    //create a VBS with fixed block dimension (see input)
    int vbmat_blocks_fmt = 1;
    int vbmat_entries_fmt = 1; //cuda needs column-major matrices
    VBS vbmat_A;

    A_row_part = linspan(0, A_rows, block_size); //row and column partitions
    A_col_part = linspan(0, A_cols, block_size);

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

    if (verbose > 1)
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

    if (verbose > 1)
    {
        cout << "VBS matrix (no zero blocks mode ON) created:" << endl;
        matprint(vbmat_A_full);
    }

    //create a VBS which is permuted with the asymmetric angle method
    VBS vbmat_A_angle;

    angle_hash_method(cmat_A, eps, A_col_part, block_cols, vbmat_A_angle, vbmat_blocks_fmt, vbmat_entries_fmt, 0);
    if (verbose > 1)
    {
        cout << "VBS matrix (asymmetric angle method) created:" << endl;
        matprint(vbmat_A_angle);
    }

    //report on block structure
    float VBS_effective_sparsity = ((float)vbmat_A_angle.nztot) / (A_rows * A_cols);
    output_couple(output_names, output_values, "VBS_AHA_effective_sparsity", VBS_effective_sparsity);
    output_couple(output_names, output_values, "VBS_AHA_block_rows", vbmat_A_angle.block_rows);

    //*******************************************
    //         MULTIPLICATION PHASE
    //___________________________________________
    //several ways of multiplying the sparse matrix
    //with a dense one, with benchmarks
    //******************************************

    //keeps track of time
    float dt;
    vec_d algo_times;
    float mean_time;
    float std_time;
   
    //output format
    cout << fixed;

    output_couple(output_names, output_values, "Warmup", warmup);
    output_couple(output_names, output_values, "Repetitions", experiment_reps);
    output_couple(output_names, output_values, "Algorithm", algo);


    if (verbose > 0)        cout << "\n \n ************************** \n STARTING THE MULTIPLICATION PHASE \n" << endl;

    //creating a random matrix X
    int B_rows = A_cols;
    int mat_B_fmt = 1;

    DataT mat_B[B_rows * B_cols] = { 0 };
    random_mat(mat_B, B_rows, B_cols, B_sparsity);

    if (verbose > 0)        std::cout << "Random matrix B created:" << std::endl;
    if (verbose > 1)        matprint(mat_B, B_rows, B_cols, B_rows, mat_B_fmt);

    //defining the output matrix C
	int C_rows = A_rows;
	int C_cols = B_cols;


    //--------------------------------------------
    //  dense-dense cublas gemm multiplication
    //--------------------------------------------
    if ((algo == 1) or (algo == -1))
    {
        //create a dense array matrix from cmat_A

        DataT mat_A[A_rows * A_cols] = { 0 };
        convert_to_mat(cmat_A, mat_A, mat_A_fmt);

        DataT mat_Cgemm[C_rows * C_cols] = { 0 };
        int mat_Cgemm_fmt = 1;

        algo_times.clear();
        for (int i = -warmup; i < experiment_reps; i++)
        {
            cublas_gemm_custom(mat_A, A_rows, A_cols, A_rows, mat_B, B_cols, B_rows, mat_Cgemm, C_rows, 1.0f, 0.0f, dt);
            //only saves non-warmup runs
            if(i >= 0) algo_times.push_back(dt);
        }

        mean_time = mean(algo_times);
        std_time = std_dev(algo_times);
        output_couple(output_names, output_values, "gemm_mean", mean_time);
        output_couple(output_names, output_values, "gemm_std", std_time);


        if (verbose > 0)        cout << "Dense-Dense multiplication. Time taken: " << mean_time << endl;
        if (verbose > 1)
        {
            cout << "GEMM Matrix:" << endl;
            matprint(mat_Cgemm, C_rows, C_cols, C_rows, 1);
        }

    }


    //--------------------------------------------
    //      VBS x dense cublas multiplication	
    //--------------------------------------------
    if ((algo == 2) or (algo == -1))
    {
        DataT mat_Cblock[C_rows * C_cols];
        int mat_Cblock_fmt = 1;

        algo_times.clear();
        for (int i = -warmup; i < experiment_reps; i++)
        {
            cublas_blockmat_multiply(vbmat_A, mat_B, B_cols, B_rows, mat_Cblock, C_rows, dt);
            //only saves non-warmup runs
            if (i >= 0) algo_times.push_back(dt);
        }

        mean_time = mean(algo_times);
        std_time = std_dev(algo_times);
        output_couple(output_names, output_values, "VBSmm_mean", mean_time);
        output_couple(output_names, output_values, "VBSmm_std", std_time);


        if (verbose > 0)
        {
            cout << "BlockSparse-Dense multiplication. Time taken: " << mean_time << endl;
        }
        if (verbose > 1)
        {

            cout << "BLOCK RESULT" << endl;
            matprint(mat_Cblock, C_rows, C_cols, C_rows, 1);
        }

        /*
        if (correct_check)
        {
            int block_success = equal(C_rows, C_cols, mat_Cgemm, C_rows, mat_Cgemm_fmt, mat_Cblock, C_rows, mat_Cblock_fmt, precision);
            if (!block_success)
            {
                std::cout << "WARNING: Block matrix multiplication test: FAILED" << std::endl;
            }
        }
        */

    }


    //--------------------------------------------
    //      VBS x dense cublas multiplication (no zero blocks mode)
    //--------------------------------------------
    if ((algo == 3) or (algo == -1))
    {
        DataT mat_Cblock_full[C_rows * C_cols];
        int mat_Cblock_full_fmt = 1;


        algo_times.clear();
        for (int i = -warmup; i < experiment_reps; i++)
        {
            cublas_blockmat_multiply(vbmat_A_full, mat_B, B_cols, B_rows, mat_Cblock_full, C_rows, dt);
            if (i >= 0) algo_times.push_back(dt);
        }

        mean_time = mean(algo_times);
        std_time = std_dev(algo_times);
        output_couple(output_names, output_values, "VBSmm_nozeros_mean", mean_time);
        output_couple(output_names, output_values, "VBSmm_nozeros_std", std_time);


        if (verbose > 0)
        {
            cout << "BlockSparse-Dense multiplication (no zero mode ON). Time taken: " << mean_time << endl;
        }
        if (verbose > 1)
        {

            cout << "BLOCK RESULT (no zero mode ON)" << endl;
            matprint(mat_Cblock_full, C_rows, C_cols, C_rows, 1);
        }

        /*
        if (correct_check)
        {
            int block_full_success = equal(C_rows, C_cols, mat_Cgemm, C_rows, mat_Cgemm_fmt, mat_Cblock_full, C_rows, mat_Cblock_full_fmt, precision);
            if (!block_full_success)
            {
                std::cout << "WARNING: Block matrix multiplication test: FAILED" << std::endl;
            }
        }
        */
    }


    //--------------------------------------------
    //      VBS x dense cublas multiplication (permuted with angle algorithm)
    //--------------------------------------------
    if ((algo == 4) or (algo == -1))
    {

        DataT mat_Cblock_angle[C_rows * C_cols];
        int mat_Cblock_angle_fmt = 1;

        algo_times.clear();
        for (int i = -warmup; i < experiment_reps; i++)
        {
            cublas_blockmat_multiply(vbmat_A_angle, mat_B, B_cols, B_rows, mat_Cblock_angle, C_rows, dt);
            if (i >= 0) algo_times.push_back(dt);
        }

        mean_time = mean(algo_times);
        std_time = std_dev(algo_times);
        output_couple(output_names, output_values, "VBSmm_mean", mean_time);
        output_couple(output_names, output_values, "VBSmm_std", std_time);

        if (verbose > 0)
        {
            cout << "BlockSparse-Dense multiplication (permuted with AHS). Time taken: " << mean_time << endl;
        }
        if (verbose > 1)
        {

            cout << "BLOCK RESULT (permuted with AHS)" << endl;
            matprint(mat_Cblock_angle, C_rows, C_cols, C_rows, mat_Cblock_angle_fmt);
        }
    }


    //--------------------------------------------
    //      CSR x Dense cusparse multiplication
    //--------------------------------------------
    if ((algo == 5) or (algo == -1))
    {

        DataT mat_C_csrmm[C_rows * C_cols];
        int mat_C_csrmm_fmt = 1;

        //prepare the cusparse CSR format
        int nnz = count_nnz(cmat_A);
        int csrRowPtr[cmat_A.rows + 1];
        int csrColInd[nnz];
        float csrVal[nnz];
        prepare_cusparse_CSR(cmat_A, csrRowPtr, csrColInd, csrVal);
        
        algo_times.clear();
        for (int i = -warmup; i < experiment_reps; i++)
        {
            cusparse_gemm_custom(cmat_A.rows, cmat_A.cols, nnz, csrRowPtr, csrColInd, csrVal, mat_B, B_cols, B_rows, mat_C_csrmm, C_rows, 1.0f, 0.0f, dt);
            if (i >= 0) algo_times.push_back(dt);
        }

        mean_time = mean(algo_times);
        std_time = std_dev(algo_times);
        output_couple(output_names, output_values, "VBSmm_mean", mean_time);
        output_couple(output_names, output_values, "VBSmm_std", std_time);


        if (verbose > 0)
        {
            cout << "CSR-Dense cusparse multiplication. Time taken: " << mean_time << endl;
        }
        if (verbose > 1)
        {

            cout << "CSR-dense cusparse:" << endl;
            matprint(mat_C_csrmm, C_rows, C_cols, C_rows, 1);
        }

        /*
        //correctness check
        int csrmm_success = equal(C_rows, C_cols, mat_Cgemm, C_rows, mat_Cgemm_fmt, mat_C_csrmm, C_rows, mat_C_csrmm_fmt, precision);
        if (!csrmm_success)
        {
            std::cout << "WARNING: CSR-Dense cusparse multiplication test: FAILED" << std::endl;
        }
        */

    }


    //cleaning
    cleanVBS(vbmat_A);
    cleanVBS(vbmat_A_full);
    cleanCSR(cmat_A);

    //OUTPUT PHASE
    if (verbose == -1)
    {
        cout << output_names << endl;
        cout << output_values << endl;
    }



}

//TODO CSR-dense cusparse multiplication

 
 