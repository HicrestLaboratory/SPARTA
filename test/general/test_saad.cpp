#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <unistd.h> //getopt, optarg
#include <ctime>
#include <cmath>
#include <algorithm> //min_element, max_element

#include <string>

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
        total_area_vec.clear();
        block_rows_vec.clear();
        nz_blocks_vec.clear();
        avg_height_vec.clear();
        skip_vec.clear();
        comparison_vec.clear();
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
    }
};

int main(int argc, char* argv[])
{

    input_parameters params;

    get_input_params(argc, argv, params);

    params.cmat_A_fmt = 0;

    CSR cmat_A;

    int mat_B_fmt = 1;
    int vbmat_blocks_fmt = 1;
    int vbmat_entries_fmt = 1; //cuda needs column-major matrices

    string output_names;
    string output_values;
    reorder_info re_info;
    Info_Collector info_collector;

    get_input_CSR(cmat_A, params);
    output_couple_parameters(params, output_names, output_values);
    cleanCSR(cmat_A);

    //*******************************************
    //	 EXPERIMENT LOOP
    //******************************************


    for (int current_repetition = 0; current_repetition < params.experiment_reps; current_repetition++)
    {

        get_input_CSR(cmat_A, params);

        if (params.verbose > 0) cout << "INPUT ACQUIRED." << endl;
        if (params.verbose > 1) matprint(cmat_A);

        //PREPARE THE SCRAMBLED AND REORDERED VBMAT
        scramble_input(cmat_A, params);

        VBS vbmat_algo;

        saad_reordering(cmat_A, vbmat_algo, params.algo_block_size, vbmat_blocks_fmt, vbmat_entries_fmt, params, re_info);

        if (params.verbose > 0)    cout << "VBS matrix (Asymmetric Angle Method) created:" << endl;
        if (params.verbose > 1)    matprint(vbmat_algo);

        info_collector.collect_info_VBS(vbmat_algo);
        info_collector.collect_info_reordering(re_info);

        cleanCSR(cmat_A);
    }



    info_collector.output_all(output_names, output_values);

    if ((params.verbose == -1) or (params.verbose > 1))
    {
        cout << output_names << endl;
        cout << output_values << endl;
    }
}