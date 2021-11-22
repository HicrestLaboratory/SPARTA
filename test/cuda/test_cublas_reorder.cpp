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
#include "reorderings.h"
#include "comp_mats.h"
#include "input.h"


using namespace std;

int main(int argc, char* argv[]) 
{

    input_parameters params;

    get_input_params(argc, argv, params);

    params.cmat_A_fmt = 1;

    CSR cmat_A;

    int mat_B_fmt = 1;
    int vbmat_blocks_fmt = 1;
    int vbmat_entries_fmt = 1; //cuda needs column-major matrices

    string output_names;
    string output_values;
    reorder_info re_info;
    Info_Collector info_collector;
    output_couple_parameters(params, output_names, output_values);

    //*******************************************
    //	 EXPERIMENT LOOP
    //******************************************


    for (int current_repetition = 0; current_repetition < params.experiment_reps; current_repetition++)
    {

        get_input_CSR(cmat_A, params);

        if (params.verbose > 0) cout << "INPUT ACQUIRED." << endl;
        if (params.verbose > 1) matprint(cmat_A);

        //PREPARE THE PERFECTLY-BLOCKED VBS
        VBS vbmat_perfect;

        convert_to_VBS(cmat_A,
            vbmat_perfect,
            params.block_size,
            params.block_size,
            vbmat_blocks_fmt, vbmat_entries_fmt);

        if (params.verbose > 0) cout << "VBS matrix created." << endl;
        if (params.verbose > 1) matprint(vbmat_perfect);


        //PREPARE THE SCRAMBLED AND REORDERED VBMAT
        scramble_input(input_cmat, params);
        
        VBS vbmat_algo;

        saad_reordering(cmat_A, vbmat_algo, params.algo_block_size, vbmat_blocks_fmt, vbmat_entries_fmt, params, re_info);

        if (params.verbose > 0)    cout << "VBS matrix (Asymmetric Angle Method) created:" << endl;
        if (params.verbose > 1)    matprint(vbmat_algo);

        info_collector.collect_info_VBS(vbmat_algo);
        info_collector.collect_info_reordering(re_info);

        //*******************************************
        //         MULTIPLICATION PHASE
        //___________________________________________
        //several ways of multiplying the sparse matrix
        //with a dense one, with benchmarks
        //******************************************

        intT A_rows = params.A_rows;
        intT A_cols = params.A_cols;
        intT B_rows = A_cols;
        intT B_cols = params.B_cols;
        intT C_rows = A_rows;
        intT C_cols = B_cols;
        int mat_B_fmt = 1;

        if (params.verbose > 0)        cout << "\n \n ************************** \n STARTING THE MULTIPLICATION PHASE \n" << endl;

        DataT* mat_B = new DataT[B_rows * B_cols]{ 0 };
        random_mat(mat_B, B_rows, B_cols, params.B_density); // creates a random DataT matrix filled with 1.000 at a fixed density

        if (params.verbose > 0)        std::cout << "Random matrix B created:" << std::endl;
        if (params.verbose > 1)        matprint(mat_B, B_rows, B_cols, B_rows, mat_B_fmt);

        //--------------------------------------------
        //      VBS perfect x dense cublas multiplication	
        //--------------------------------------------

        if (params.verbose > 0)        cout << "Starting VBS-dense cublas multiplication" << endl;

        DataT* mat_Cperfect_block = new DataT[C_rows * C_cols];

        for (int i = -params.warmup; i < 1; i++)//do warmup runs
        {
            float dt = 0;
            cublas_blockmat_multiply(vbmat_algo, mat_B, B_cols, B_rows, mat_Cperfect_block, C_rows, dt, params.n_streams);
            if (i >= 0) info_collector.vbs_perfect_times.push_back(dt);
            if (params.verbose > 0)            cout << "BlockSparse-Dense multiplication. Time taken(ms): " << dt << endl;
        }

        delete[] mat_Cperfect_block;
        cleanVBS(vbmat_algo);



        //--------------------------------------------
        //      VBS algo x dense cublas multiplication	
        //--------------------------------------------

        if (params.verbose > 0)        cout << "Starting VBS-dense cublas multiplication" << endl;

        DataT* mat_Cblock = new DataT[C_rows * C_cols];

        for (int i = -params.warmup; i < 1; i++)//do warmup runs
        {
            float dt = 0;
            cublas_blockmat_multiply(vbmat_perfect, mat_B, B_cols, B_rows, mat_Cblock, C_rows, dt, params.n_streams);
            if (i >= 0) info_collector.vbs_algo_times.push_back(dt);
            if (params.verbose > 0)            cout << "BlockSparse-Dense multiplication. Time taken(ms): " << dt << endl;
        }

        delete[] mat_Cblock;
        cleanVBS(vbmat_perfect);

        //--------------------------------------------
        //      CSR x Dense cusparse multiplication
        //--------------------------------------------
        if (typeid(DataT) != typeid(float)) cout << "WARNING: only float supported for CUSPARSE. DataT can be changed in sparse_utilities.h" << endl;

        if (params.verbose > 0) cout << "Starting cusparse-dense cublas multiplication" << endl;

        DataT* mat_C_csrmm = new DataT[C_rows * C_cols];
        int mat_C_csrmm_fmt = 1;

        //prepare the cusparse CSR format
        int A_nnz = params.A_nnz;
        int* csrRowPtr = new int[A_rows + 1];
        int* csrColInd = new int[A_nnz];
        float* csrVal = new float[A_nnz];
        prepare_cusparse_CSR(cmat_A, csrRowPtr, csrColInd, csrVal);

        for (int i = -params.warmup; i < 1; i++)
        {
            float dt = 0;
            cusparse_gemm_custom(A_rows, A_cols, A_nnz, csrRowPtr, csrColInd, csrVal, mat_B, B_cols, B_rows, mat_C_csrmm, C_rows, 1.0f, 0.0f, dt);
            if (i >= 0) info_collector.cusparse_times.push_back(dt);
            if (params.verbose > 0)          cout << "CSR-Dense cusparse multiplication. Time taken: " << dt << endl;
        }

        delete[] mat_C_csrmm;
        delete[] csrColInd;
        delete[] csrRowPtr;
        delete[] csrVal;
        cleanCSR(cmat_A);

        delete[] mat_B;
    }




    //*******************************************
    //	 EXPERIMENT LOOP
    //******************************************

    info_collector.output_all(output_names, output_values);

    if ((params.verbose == -1) or (params.verbose > 1))
    {
        cout << output_names << endl;
        cout << output_values << endl;
    }
}
 
 