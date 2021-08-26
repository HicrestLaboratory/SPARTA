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
#include "input.h"


using namespace std;

int main(int argc, char* argv[]) {

    input_parameters params;

    get_input_params(argc, argv, params);

    params.cmat_A_fmt = 1;

    CSR input_cmat; //this will hold the CSR matrix

    get_input_CSR(input_cmat, params);

    srand(params.seed);
    //_____________________________
    //*******************************************
    //		END OF INPUT
    //cmat must hold a proper CSR matrix at this point
    //******************************************

    if (params.verbose > 0) cout << "INPUT ACQUIRED." << endl;
    if (params.verbose > 1) matprint(input_cmat);

    //update rows and cols count to input values

    string output_names;
    string output_values;

    output_couple_parameters(params, output_names, output_values);

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
    if (params.verbose > 0)        cout << "\n \n ************************** \n STARTING THE MULTIPLICATION PHASE \n" << endl;

    //creating a random matrix X
    intT B_rows = params.A_cols;
    int mat_B_fmt = 1;

    //TODO smart pointers for matrices

    DataT* mat_B = new DataT[B_rows * params.B_cols]{ 0 };
    random_mat(mat_B, B_rows, params.B_cols, params.B_density); // creates a random DataT matrix filled with 1.000 at a fixed density

    if (params.verbose > 0)        cout << "Random matrix B created:" << endl;
    if (params.verbose > 1)        matprint(mat_B, B_rows, params.B_cols, B_rows, mat_B_fmt);

    //defining the output matrix C
	int C_rows = params.A_rows;
	int C_cols = params.B_cols;

    //--------------------------------------------
    //  dense-dense cublas gemm multiplication
    //--------------------------------------------
    DataT* mat_A_gemm = new DataT [params.A_rows * params.A_cols]{ 0 };
 
    convert_to_mat(input_cmat, mat_A_gemm, params.cmat_A_fmt);
 
    DataT* mat_Cgemm = new DataT[C_rows * C_cols]{ 0 };
    int mat_Cgemm_fmt = 1;

    algo_times.clear();
 
    for (int i = -params.warmup; i < params.experiment_reps; i++)
    {
        cublas_gemm_custom(mat_A_gemm, params.A_rows, params.A_cols, params.A_rows, mat_B, params.B_cols, B_rows, mat_Cgemm, C_rows, 1.0f, 0.0f, dt);
        //only saves non-warmup runs
        if(i >= 0) algo_times.push_back(dt);
    }

    mean_time = mean(algo_times);
    std_time = std_dev(algo_times);
    output_couple(output_names, output_values, "gemm_mean(ms)", mean_time);
    output_couple(output_names, output_values, "gemm_std", std_time);

    if (params.verbose > 0)        cout << "Dense-Dense multiplication. Time taken(ms): " << mean_time << endl;
    if (params.verbose > 1)
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
        int* csrRowPtr = new int[input_cmat.rows + 1];
        int* csrColInd = new int[params.A_nnz];
        float* csrVal = new float[params.A_nnz];
        prepare_cusparse_CSR(input_cmat, csrRowPtr, csrColInd, csrVal);
        
        algo_times.clear();
        for (int i = -params.warmup; i < params.experiment_reps; i++)
        {
            cusparse_gemm_custom(params.A_rows, params.A_cols, params.A_nnz, csrRowPtr, csrColInd, csrVal, mat_B, params.B_cols, B_rows, mat_C_csrmm, C_rows, 1.0f, 0.0f, dt);
            if (i >= 0) algo_times.push_back(dt);
        }

        mean_time = mean(algo_times);
        std_time = std_dev(algo_times);
        output_couple(output_names, output_values, "cusparse_spmm_mean(ms)", mean_time);
        output_couple(output_names, output_values, "cusparse_spmm_std", std_time);


        if (params.verbose > 0)
        {
            cout << "CSR-Dense cusparse multiplication. Time taken: " << mean_time << endl;
        }
        if (params.verbose > 1)
        {

            cout << "CSR-dense cusparse:" << endl;
            matprint(mat_C_csrmm, C_rows, C_cols, C_rows, 1);
        }

        delete[] mat_C_csrmm;
        delete[] csrColInd;
        delete[] csrRowPtr;
        delete[] csrVal;

    }


//--------------------------------------------
//      OUTPUT PHASE
//--------------------------------------------
    if ((params.verbose == -1) or (params.verbose > 1))
    {
        cout << output_names << endl;
        cout << output_values << endl;
    }

//--------------------------------------------
//      CLEANING
//--------------------------------------------
    delete[] mat_B;
    cleanCSR(input_cmat);
}
 
 
