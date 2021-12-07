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

struct Info_Collector
{
    svi total_area_vec;
    svi block_rows_vec;
    svi nz_blocks_vec;
    vec_d avg_height_vec;
    vec_d skip_vec;
    vec_d comparison_vec;
    vec_d vbs_algo_times;
    vec_d vbs_perfect_times;
    vec_d cusparse_times;
    string output_names = "";
    string output_values = "";

    void collect_info_VBS(VBS& vbmat)
    {
        //collect info for a post-reordering VBMAT

        //accumulate results in vectors
        avg_height_vec.push_back(vbmat.get_avg_height());
        total_area_vec.push_back(vbmat.nztot);
        block_rows_vec.push_back(vbmat.block_rows);
        nz_blocks_vec.push_back(vbmat.nz_blocks());
    }

    void collect_info_reordering(reorder_info re_info)
    {
        skip_vec.push_back(re_info.skipped);
        comparison_vec.push_back(re_info.comparisons);
    }

    void clean()
    {
        output_names = "";
        output_values = "";
        total_area_vec.clear();
        block_rows_vec.clear();
        nz_blocks_vec.clear();
        avg_height_vec.clear();
        skip_vec.clear();
        comparison_vec.clear();
        vbs_algo_times.clear();
        cusparse_times.clear();
        vbs_perfect_times.clear();
    }

    void output_all(string& output_names, string& output_values)
    {
        output_couple(output_names, output_values, "VBS_avg_nzblock_height", mean(avg_height_vec));
        output_couple(output_names, output_values, "VBS_avg_nzblock_height_error", std_dev(avg_height_vec));

        output_couple(output_names, output_values, "VBS_total_nonzeros", mean(total_area_vec));
        output_couple(output_names, output_values, "VBS_total_nonzeros_error", std_dev(total_area_vec));

        output_couple(output_names, output_values, "VBS_block_rows", mean(block_rows_vec));
        output_couple(output_names, output_values, "VBS_block_rows_error", std_dev(block_rows_vec));

        output_couple(output_names, output_values, "VBS_nz_blocks", mean(nz_blocks_vec));
        output_couple(output_names, output_values, "VBS_nz_blocks_error", std_dev(nz_blocks_vec));

        output_couple(output_names, output_values, "avg_skipped", mean(skip_vec));
        output_couple(output_names, output_values, "skipped_std", std_dev(skip_vec));

        output_couple(output_names, output_values, "avg_comparisons", mean(comparison_vec));
        output_couple(output_names, output_values, "comparisons_std", std_dev(comparison_vec));

        output_couple(output_names, output_values, "VBSmm_algo_mean(ms)", mean(vbs_algo_times));
        output_couple(output_names, output_values, "VBSmm_algo_std", std_dev(vbs_algo_times));

        if (!vbs_perfect_times.empty())
        {
            output_couple(output_names, output_values, "VBSmm_perfect_mean(ms)", mean(vbs_perfect_times));
            output_couple(output_names, output_values, "VBSmm_perfect_std", std_dev(vbs_perfect_times));
        }

        output_couple(output_names, output_values, "cusparse_spmm_mean(ms)", mean(cusparse_times));
        output_couple(output_names, output_values, "cusparse_spmm_std", std_dev(cusparse_times));

    }
};

int main(int argc, char* argv[]) 
{

    input_parameters params;

    get_input_params(argc, argv, params);


    params.exp_name = "reorder_and_multiplication_optimized";
    params.cmat_A_fmt = 1;

    CSR cmat_A;

    int mat_B_fmt = 1;
    int vbmat_blocks_fmt = 1;
    int vbmat_entries_fmt = 1; //cuda needs column-major matrices

    string output_names;
    string output_values;
    reorder_info re_info;
    Info_Collector info_collector;

    get_input_CSR(cmat_A, params);

    if (params.verbose > 0) cout << "INPUT ACQUIRED." << endl;
    if (params.verbose > 1) matprint(cmat_A);


    //*******************************************
    //	 EXPERIMENT LOOP
    //******************************************


    vector<intT> Ns = { 4096, 8192, 16384, 32768 };

    //PREPARE THE PERFECTLY-BLOCKED VBS
    VBS vbmat_perfect;
        
    if (params.algo == 1)
    {
        convert_to_VBS(cmat_A,
            vbmat_perfect,
            params.block_size,
            params.block_size,
            vbmat_blocks_fmt, vbmat_entries_fmt);

        if (params.verbose > 0) cout << "VBS perfect matrix created." << endl;
        if (params.verbose > 1) matprint(vbmat_perfect);
    }

    //PREPARE THE SCRAMBLED AND REORDERED VBMAT
    scramble_input(cmat_A, params);
        
    VBS vbmat_algo;

    saad_reordering(cmat_A, vbmat_algo, params.algo_block_size, vbmat_blocks_fmt, vbmat_entries_fmt, params, re_info);

    if (params.verbose > 0)    cout << "VBS matrix (Asymmetric Angle Method) created:" << endl;
    if (params.verbose > 1)    matprint(vbmat_algo);

    for (intT N : Ns)
    {
        params.B_cols = N;
        output_couple_parameters(params, info_collector.output_names, info_collector.output_values);
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


        if (params.algo == 1)
        {
            if (params.verbose > 0)        cout << "Starting VBS-perfect-dense cublas multiplication" << endl;

            DataT_C* mat_Cblock = new DataT_C[C_rows * C_cols];

            for (int i = -params.warmup; i < params.experiment_reps; i++)//do warmup runs
            {
                float dt = 0;
                cublas_blockmat_multiply(vbmat_perfect, mat_B, B_cols, B_rows, mat_Cblock, C_rows, dt, params.n_streams);
                if (i >= 0) info_collector.vbs_perfect_times.push_back(dt);
                if (params.verbose > 0)            cout << "BlockSparse-Dense multiplication. Time taken(ms): " << dt << endl;
            }

            delete[] mat_Cblock;
        }


        //--------------------------------------------
        //      VBS algo x dense cublas multiplication	
        //--------------------------------------------

        if (params.verbose > 0)        cout << "Starting VBS-reordered-dense cublas multiplication" << endl;

        DataT_C* mat_Cblock = new DataT_C[C_rows * C_cols];

        for (int i = -params.warmup; i < params.experiment_reps; i++)//do warmup runs
        {
            float dt = 0;
            cublas_blockmat_multiply(vbmat_algo, mat_B, B_cols, B_rows, mat_Cblock, C_rows, dt, params.n_streams);
            if (i >= 0) info_collector.vbs_algo_times.push_back(dt);
            if (params.verbose > 0)            cout << "BlockSparse-Dense multiplication. Time taken(ms): " << dt << endl;
        }

        delete[] mat_Cblock;

        //--------------------------------------------
        //      CSR x Dense cusparse multiplication
        //--------------------------------------------

        if (params.verbose > 0) cout << "Starting cusparse-dense cublas multiplication" << endl;

        DataT_C* mat_C_csrmm = new DataT_C[C_rows * C_cols];
        int mat_C_csrmm_fmt = 1;

        //prepare the cusparse CSR format
        int A_nnz = params.A_nnz;
        int* csrRowPtr = new int[A_rows + 1];
        int* csrColInd = new int[A_nnz];
        DataT* csrVal = new DataT[A_nnz];
        prepare_cusparse_CSR(cmat_A, csrRowPtr, csrColInd, csrVal);

        for (int i = -params.warmup; i < params.experiment_reps; i++)
        {
            float dt = 0;
            cusparse_gemm_custom(A_rows, A_cols, A_nnz, csrRowPtr, csrColInd, csrVal, mat_B, B_cols, B_rows, mat_C_csrmm, C_rows, 1, 0, dt);
            if (i >= 0) info_collector.cusparse_times.push_back(dt);
            if (params.verbose > 0)          cout << "CSR-Dense cusparse multiplication. Time taken: " << dt << endl;
        }

        delete[] mat_C_csrmm;
        delete[] csrColInd;
        delete[] csrRowPtr;
        delete[] csrVal;

        delete[] mat_B;




        //*******************************************
        //	 EXPERIMENT LOOP
        //******************************************

        info_collector.output_all(output_names, output_values);

        if ((params.verbose == -1) or (params.verbose > 1))
        {
            cout << info_collector.output_names << endl;
            cout << info_collector.output_values << endl;
        }

        info_collector.clean();
    }


    if (params.algo == 1) cleanVBS(vbmat_perfect);
    cleanVBS(vbmat_algo);
    cleanCSR(cmat_A);
}
 
 