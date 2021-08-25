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
    //mean of a vector
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
    //std of a vector
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
    //append name to names and value to values; add spaces
    //used to produce output in CSV-like form;

    names += name + " ";
    values += to_string(value) + " ";
}

int main(int argc, char* argv[]) {

    opterr = 0;

    int verbose = 3;

    int input_type = 1;
    int A_rows = 12;            //rows in the square input matrix;
    int A_cols = 8;
    int mat_A_fmt = 1;          //cuda needs column-major matrices
    float density = 0.5;        //density of the input matrix;
    string input_source;

    int B_cols = 5;             //number of columns in the output matrix;
    float B_density = 1.;       //density of the multiplication matrix

    int seed = 123;

    int warmup = 1;             //number of warmup experiments
    int experiment_reps = 5;    //number of non-warmup repetitions
    bool header = 1;
	
    //terminal options loop
    opterr = 0;
    char c;
    while ((c = getopt(argc, argv, "i:q:f:m:n:r:k:v:w:S:H")) != -1)
        switch (c)
        {
        case 'i':// select input example
            input_type = stoi(optarg);
            //  1: Random CSR
            //  2: SNAP Edgelist
            break;

        case 'f': //select source file
            //has only effect for example 2 and 3;
            input_source = optarg;
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

        case 'q': //density of input matrix
            //has only effect for example 1 and 4
            density = stof(optarg);
            if (density < 0 or density > 1) {
                fprintf(stderr, "Option -k tried to set density outside of [0,1]");
                return 1;
            }
            break;

        case 'r': //number of experiment repetitions
            experiment_reps = stoi(optarg);
            break;

	case 'S': //random seed
            seed = stoi(optarg);
            break;

        case 'v': //verbose
            verbose = stoi(optarg);
            break;

        case 'w': //warmup repetitions
            warmup = stoi(optarg);
            break;
	case 'H':
	    header = stoi(optarg);
	    break;
        case '?':
            fprintf(stderr, "Option -%c does not exists, or requires an argument.\n", optopt);
            return 1;
        default:
            abort();
        }

    //TODO -h HELP

    srand(seed);

    //INPUT CONVERSION TO Compressed Sparse Row (CSR)

    CSR cmat_A; //this will hold the CSR matrix
    int cmat_A_fmt = 0;


    //INPUT EXAMPLE 1: RANDOM CSR
    //create a random sparse matrix
    if (input_type == 1) {
        DataT* rand_mat = new DataT[A_cols * A_rows];
        random_mat(rand_mat, A_rows, A_cols, density); //generate random mat

        convert_to_CSR(rand_mat, A_rows, A_cols, mat_A_fmt, cmat_A, cmat_A_fmt);
        delete[] rand_mat;

        if (verbose > 0) cout << "CREATED A RANDOM CSR with density = " << density << endl;
    }
    //______________________________________


    //TEST
    //INPUT EXAMPLE 2: read graph in edgelist format into CSR
    else if (input_type == 2)
    {
        //TEST
        //INPUT EXAMPLE 2: read graph in edgelist format into CSR
        if (input_source.empty()) input_source = "testgraph1.txt";

        string delimiter = "\t";
        read_edgelist(input_source, input_cmat, input_cmat_fmt, delimiter);
        if (verbose > 0) cout << "IMPORTED A CSR FROM A SNAP EDGELIST" << endl;
        //______________________________________
    }
     //______________________________________

    //___________________________________________
    //*******************************************
    //		END OF INPUT
    //cmat must hold a proper CSR matrix at this point
    //******************************************

    if (verbose > 0) cout << "INPUT ACQUIRED." << endl;
    if (verbose > 1) matprint(cmat_A);

    //update rows and cols count to input values
    A_rows = cmat_A.rows;
    A_cols = cmat_A.cols;

    string output_names;
    string output_values;

    int A_nnz = count_nnz(cmat_A);

    output_couple(output_names, output_values, "input_type", input_type);
    output_couple(output_names, output_values, "A_rows", cmat_A.rows);
    output_couple(output_names, output_values, "A_cols", cmat_A.cols);
    output_couple(output_names, output_values, "A_total_nonzeros", A_nnz);
    output_couple(output_names, output_values, "A_entries_density", density);

    output_couple(output_names, output_values, "B_cols", B_cols);
    output_couple(output_names, output_values, "B_density", B_density);


    //******************************************
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


    if (verbose > 0)        cout << "\n \n ************************** \n STARTING THE MULTIPLICATION PHASE \n" << endl;

    //creating a random matrix X
    int B_rows = A_cols;
    int mat_B_fmt = 1;


    //TODO smart pointers for matrices

    DataT* mat_B = new DataT[B_rows * B_cols]{ 0 };
    random_mat(mat_B, B_rows, B_cols, B_density); // creates a random DataT matrix filled with 1.000 at a fixed density

    if (verbose > 0)        std::cout << "Random matrix B created:" << std::endl;
    if (verbose > 1)        matprint(mat_B, B_rows, B_cols, B_rows, mat_B_fmt);

    //defining the output matrix C
	int C_rows = A_rows;
	int C_cols = B_cols;

    //--------------------------------------------
    //  dense-dense cublas gemm multiplication
    //--------------------------------------------
        DataT* mat_A_gemm = new DataT [A_rows * A_cols]{ 0 };
 
        convert_to_mat(cmat_A, mat_A_gemm, mat_A_fmt);
 
        DataT* mat_Cgemm = new DataT[C_rows * C_cols]{ 0 };
        int mat_Cgemm_fmt = 1;

        algo_times.clear();
 
        for (int i = -warmup; i < experiment_reps; i++)
        {
            cublas_gemm_custom(mat_A_gemm, A_rows, A_cols, A_rows, mat_B, B_cols, B_rows, mat_Cgemm, C_rows, 1.0f, 0.0f, dt);
            //only saves non-warmup runs
            if(i >= 0) algo_times.push_back(dt);
        }

        mean_time = mean(algo_times);
        std_time = std_dev(algo_times);
        output_couple(output_names, output_values, "gemm_mean(ms)", mean_time);
        output_couple(output_names, output_values, "gemm_std", std_time);

        if (verbose > 0)        cout << "Dense-Dense multiplication. Time taken(ms): " << mean_time << endl;
        if (verbose > 1)
        {
            cout << "GEMM Matrix:" << endl;
            matprint(mat_Cgemm, C_rows, C_cols, C_rows, 1);
        }
        
        delete[] mat_A_gemm;
        delete[] mat_Cgemm;

    //--------------------------------------------
    //      CSR x Dense cusparse multiplication
    //--------------------------------------------
        if (typeid(DataT) != typeid(float))
        {
            if (verbose > 0)         cout << "WARNING: only float supported for CUSPARSE. DataT can be changed in sparse_utilities.h" << endl;
        }
        else
        {

            DataT* mat_C_csrmm = new DataT[C_rows * C_cols];
            int mat_C_csrmm_fmt = 1;

            //prepare the cusparse CSR format
            int nnz = count_nnz(cmat_A);
            int* csrRowPtr = new int[cmat_A.rows + 1];
            int* csrColInd = new int[nnz];
            float* csrVal = new float[nnz];
            prepare_cusparse_CSR(cmat_A, csrRowPtr, csrColInd, csrVal);
        
            algo_times.clear();
            for (int i = -warmup; i < experiment_reps; i++)
            {
                cusparse_gemm_custom(cmat_A.rows, cmat_A.cols, nnz, csrRowPtr, csrColInd, csrVal, mat_B, B_cols, B_rows, mat_C_csrmm, C_rows, 1.0f, 0.0f, dt);
                if (i >= 0) algo_times.push_back(dt);
            }

            mean_time = mean(algo_times);
            std_time = std_dev(algo_times);
            output_couple(output_names, output_values, "cusparse_spmm_mean(ms)", mean_time);
            output_couple(output_names, output_values, "cusparse_spmm_std", std_time);


            if (verbose > 0)
            {
                cout << "CSR-Dense cusparse multiplication. Time taken: " << mean_time << endl;
            }
            if (verbose > 1)
            {

                cout << "CSR-dense cusparse:" << endl;
                matprint(mat_C_csrmm, C_rows, C_cols, C_rows, 1);
            }

            delete[] mat_C_csrmm;
            delete[] csrColInd;
            delete[] csrRowPtr;
            delete[] csrVal;

        }


    //cleaning

    //OUTPUT PHASE
    if ((verbose == -1) or (verbose > 1))
    {
        if(header) cout << output_names << endl;
        cout << output_values << endl;
    }

    delete[] mat_B;
    cleanCSR(cmat_A);
}
 
 
