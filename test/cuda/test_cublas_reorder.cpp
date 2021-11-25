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
    svi min_block_vec;
    svi max_block_vec;
    vec_d avg_height_vec;
    vec_d skip_vec;
    vec_d comparison_vec;
    vec_d vbs_algo_times;
    vec_d vbs_perfect_times;
    vec_d cusparse_times;

    void collect_info_VBS(VBS& vbmat)
    {
        //collect info for a post-reordering VBS

        intT max_block_H = 0;
        intT min_block_H = vbmat.rows();
        float avg_block_height = 0.;
        intT tot_nz_blocks = 0;
        for (intT i = 0; i < vbmat.block_rows; i++)
        {
            intT b_size = vbmat.row_part[i + 1] - vbmat.row_part[i];
            cout < b_size < endl;
            avg_block_height += b_size * vbmat.nzcount[i];
            tot_nz_blocks += vbmat.nzcount[i];
            if (b_size > max_block_H) max_block_H = b_size;
            if (b_size < min_block_H) min_block_H = b_size;
        }
        avg_block_height /= tot_nz_blocks;

        //accumulate results in vectors
        avg_height_vec.push_back(avg_block_height);
        total_area_vec.push_back(vbmat.nztot);
        block_rows_vec.push_back(vbmat.block_rows);
        nz_blocks_vec.push_back(tot_nz_blocks);
        min_block_vec.push_back(min_block_H);
        max_block_vec.push_back(max_block_H);
    }

    void collect_info_reordering(reorder_info re_info)
    {
        skip_vec.push_back(re_info.skipped);
        comparison_vec.push_back(re_info.comparisons);
    }

    void clean()
    {
        total_area_vec.clear();
        block_rows_vec.clear();
        nz_blocks_vec.clear();
        min_block_vec.clear();
        max_block_vec.clear();
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

        output_couple(output_names, output_values, "VBS_min_block_H", mean(min_block_vec));
        output_couple(output_names, output_values, "VBS_min_block_H_error", std_dev(min_block_vec));

        output_couple(output_names, output_values, "VBS_max_block_H", mean(max_block_vec));
        output_couple(output_names, output_values, "VBS_max_block_H_error", std_dev(max_block_vec));

        output_couple(output_names, output_values, "avg_skipped", mean(skip_vec));
        output_couple(output_names, output_values, "skipped_std", std_dev(skip_vec));

        output_couple(output_names, output_values, "avg_comparisons", mean(comparison_vec));
        output_couple(output_names, output_values, "comparisons_std", std_dev(comparison_vec));


        output_couple(output_names, output_values, "VBSmm_algo_mean(ms)", mean(vbs_algo_times));
        output_couple(output_names, output_values, "VBSmm_algo_std", std_dev(vbs_algo_times));

        output_couple(output_names, output_values, "VBSmm_perfect_mean(ms)", mean(vbs_perfect_times));
        output_couple(output_names, output_values, "VBSmm_perfect_std", std_dev(vbs_perfect_times));

        output_couple(output_names, output_values, "cusparse_spmm_mean(ms)", mean(cusparse_times));
        output_couple(output_names, output_values, "cusparse_spmm_std", std_dev(cusparse_times));

    }
};

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

        if (params.verbose > 0) cout << "VBS perfect matrix created." << endl;
        if (params.verbose > 1) matprint(vbmat_perfect);


        //PREPARE THE SCRAMBLED AND REORDERED VBMAT
        scramble_input(cmat_A, params);
        
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
 
 