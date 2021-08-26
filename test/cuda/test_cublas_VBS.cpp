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

    CSR cmat_A; //this will hold the CSR matrix

    get_input_CSR(cmat_A, params);

    srand(params.seed);

    //___________________________________________
    //*******************************************
    //		END OF INPUT
    //cmat must hold a proper CSR matrix at this point
    //******************************************

    if (params.verbose > 0) cout << "INPUT ACQUIRED." << endl;
    if (params.verbose > 1) matprint(cmat_A);

    //update rows and cols count to input values
    
    intT A_rows = params.A_rows;
    intT A_cols = params.A_cols;
    intT B_rows = A_cols;
    intT B_cols = params.B_cols;
    int mat_B_fmt = 1;

    intT C_rows = A_rows;
    intT C_cols = B_cols;

    //*******************************************
    //	 OUTPUT PARAMETERS
    //******************************************


    string output_names;
    string output_values;

    output_couple_parameters(params, output_names, output_values);

    //*******************************************
    //	 CREATE A VBS
    //******************************************

    int vbmat_blocks_fmt = 1;
    int vbmat_entries_fmt = 1; //cuda needs column-major matrices
    intT block_rows = A_rows / params.block_size;
    intT block_cols = A_cols / params.block_size;

    intT* A_row_part = new intT[block_rows + 1]; //partitions have one element more for the rightmost border.
    intT* A_col_part = new intT[block_cols + 1];
    partition(A_row_part, 0, cmat_A.rows, params.block_size); //row and column partitions (TODO make it work when block_size does not divide rows)
    partition(A_col_part, 0, cmat_A.cols, params.block_size);

    //Create a VBS with fixed block dimension (see input)
    VBS vbmat_A;
    if (params.algo == 2 or params.algo == -1)
    {

        if ((A_rows % params.block_size != 0) or (A_cols % params.block_size != 0))
        {
            std::cout << "WARNING: The row or column dimension of the input matrix is not multiple of the block size " << std::endl;
        }

        convert_to_VBS(cmat_A,
            vbmat_A,
            block_rows, A_row_part,
            block_cols, A_col_part,
            vbmat_blocks_fmt, vbmat_entries_fmt);

        if (params.verbose > 0) cout << "VBS matrix created." << endl;
        if (params.verbose > 1) matprint(vbmat_A);
    }

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


    if (params.verbose > 0)        cout << "\n \n ************************** \n STARTING THE MULTIPLICATION PHASE \n" << endl;

    //TODO smart pointers for matrices

    DataT* mat_B = new DataT[B_rows * B_cols]{ 0 };
    random_mat(mat_B, B_rows, B_cols, params.B_density); // creates a random DataT matrix filled with 1.000 at a fixed density

    if (params.verbose > 0)        std::cout << "Random matrix B created:" << std::endl;
    if (params.verbose > 1)        matprint(mat_B, B_rows, B_cols, B_rows, mat_B_fmt);

    //defining the output matrix C

    DataT* mat_Cgemm;


    //--------------------------------------------
    //  dense-dense cublas gemm multiplication
    //--------------------------------------------
    if ((params.algo == 1) or (params.algo == -1))
    {
        if (A_rows * A_cols < 10000000) //avoid this step if the matrix is too big
        {

            if (params.verbose > 0)        cout << "Starting dense-dense cublas multiplication" << endl;

            //create a dense array matrix from cmat_A

            DataT* mat_A_gemm = new DataT[A_rows * A_cols]{ 0 };

            int mat_A_fmt = 1;
            convert_to_mat(cmat_A, mat_A_gemm, mat_A_fmt);

            mat_Cgemm = new DataT[C_rows * C_cols]{ 0 };
            int mat_Cgemm_fmt = 1;

            algo_times.clear();


            for (int i = -params.warmup; i < params.experiment_reps; i++)
            {
                cublas_gemm_custom(mat_A_gemm, A_rows, A_cols, A_rows, mat_B, B_cols, B_rows, mat_Cgemm, C_rows, 1.0f, 0.0f, dt);
                //only saves non-warmup runs
                if (i >= 0) algo_times.push_back(dt);
            }

            mean_time = mean(algo_times);
            std_time = std_dev(algo_times);
            output_couple(output_names, output_values, "gemm_mean(ms)", mean_time);
            output_couple(output_names, output_values, "gemm_std", std_time);

            if (params.verbose > 0)        cout << "Dense-Dense multiplication. Time taken(ms): " << mean_time << endl;
            if (params.verbose > 1)
            {
                matprint(mat_Cgemm, C_rows, C_cols, C_rows, 1);
            }

            delete[] mat_A_gemm;
        }
    }

    //--------------------------------------------
    //      VBS x dense cublas multiplication	
    //--------------------------------------------
    if ((params.algo == 2) or (params.algo == -1))
    {

        if (params.verbose > 0)        cout << "Starting VBS-dense cublas multiplication" << endl;

        DataT* mat_Cblock = new DataT[C_rows * C_cols];
        int mat_Cblock_fmt = 1;

        algo_times.clear();
        for (int i = -params.warmup; i < params.experiment_reps; i++)
        {
            cublas_blockmat_multiply(vbmat_A, mat_B, B_cols, B_rows, mat_Cblock, C_rows, dt, params.n_streams);
            //only saves non-warmup runs
            if (i >= 0) algo_times.push_back(dt);
        }

        mean_time = mean(algo_times);
        std_time = std_dev(algo_times);
        output_couple(output_names, output_values, "VBSmm_mean(ms)", mean_time);
        output_couple(output_names, output_values, "VBSmm_std", std_time);

        if (params.check_correct)
        {
            bool vbs_check = equal(C_rows, C_cols, mat_Cgemm, C_cols, 0, mat_Cblock, C_cols, 0, 0.00001f);
            output_couple(output_names, output_values, "vbs_check", vbs_check);
        }

        if (params.verbose > 0)
        {
            cout << "BlockSparse-Dense multiplication. Time taken(ms): " << mean_time << endl;
        }
        if (params.verbose > 1)
        {

            cout << "BLOCK RESULT" << endl;
            matprint(mat_Cblock, C_rows, C_cols, C_rows, 1);
        }

        delete[] mat_Cblock;

    }

    //--------------------------------------------
    //      CSR x Dense cusparse multiplication
    //--------------------------------------------
    if ((params.algo == 5) or (params.algo == -1))
        if (typeid(DataT) != typeid(float))
        {
            if (params.verbose > 0)         cout << "WARNING: only float supported for CUSPARSE. DataT can be changed in sparse_utilities.h" << endl;
        }
        else
        {

            if (params.verbose > 0)        cout << "Starting cusparse-dense cublas multiplication" << endl;

            DataT* mat_C_csrmm = new DataT[C_rows * C_cols];
            int mat_C_csrmm_fmt = 1;

            //prepare the cusparse CSR format
            int A_nnz = params.A_nnz;
            int* csrRowPtr = new int[A_rows + 1];
            int* csrColInd = new int[A_nnz];
            float* csrVal = new float[A_nnz];
            prepare_cusparse_CSR(cmat_A, csrRowPtr, csrColInd, csrVal);
        
            algo_times.clear();
            for (int i = -params.warmup; i < params.experiment_reps; i++)
            {
                cusparse_gemm_custom(A_rows, A_cols, A_nnz, csrRowPtr, csrColInd, csrVal, mat_B, B_cols, B_rows, mat_C_csrmm, C_rows, 1.0f, 0.0f, dt);
                if (i >= 0) algo_times.push_back(dt);
            }

            mean_time = mean(algo_times);
            std_time = std_dev(algo_times);
            output_couple(output_names, output_values, "cusparse_spmm_mean(ms)", mean_time);
            output_couple(output_names, output_values, "cusparse_spmm_std", std_time);

            if (params.check_correct)
            {
                bool csr_check = equal(C_rows, C_cols, mat_Cgemm, C_cols, 0, mat_C_csrmm, C_cols, 0, 0.00001f);
                output_couple(output_names, output_values, "csr_check", csr_check);
            }

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


    //cleaning

    //OUTPUT PHASE
    if ((params.verbose == -1) or (params.verbose > 1))
    {
        cout << output_names << endl;
        cout << output_values << endl;
    }

    delete[] mat_B;
    cleanCSR(cmat_A);
    if (params.algo == 2 or params.algo == -1) cleanVBS(vbmat_A);

}
 
 