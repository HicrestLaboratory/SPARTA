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


int main(int argc, char* argv[]) {

    input_parameters params;

    get_input_params(argc, argv, params);

    CSR input_cmat; //this will hold the CSR matrix
    
    srand(params.seed);

    get_input_CSR(input_cmat, params);

    //*******************************************
    //		END OF INPUT
    //spmat must hold a proper CSR matrix at this point
    //******************************************

    if (params.verbose > 0) cout << "INPUT ACQUIRED." << endl;
    if (params.verbose > 1) matprint(input_cmat);
  
    //*******************************************
    //	 STORE PARAMETERS 
    //******************************************

    string output_names;
    string output_values;

    output_couple_parameters(params, output_names, output_values);

    //*******************************************
    //	 RUN REORDERING EXPERIMENTS 
    //******************************************

    svi total_area_vec;
    svi block_rows_vec;
    svi nz_blocks_vec;
    svi min_block_vec;
    svi max_block_vec;
    vec_d avg_height_vec;

    for (int current_repetition = 0; current_repetition < params.experiment_reps; current_repetition++)
    {
        scramble_input(input_cmat, params);

        int vbmat_blocks_fmt = 0;
        int vbmat_entries_fmt = 0;
        intT algo_block_cols = std::ceil((float)params.A_cols / params.algo_block_size);

        //prepare the column partition
        intT* algo_col_part = new intT[algo_block_cols + 1]; 
        partition(algo_col_part, 0, params.A_cols, params.algo_block_size);

        //run the reordering algo
        intT* hash_groups = new intT[params.A_rows];
        saad_reordering(input_cmat, params, hash_groups);

        //create the block matrix
        VBS vbmat_algo;
        group_to_VBS(input_cmat, hash_groups, algo_col_part, algo_block_cols, vbmat_algo, vbmat_blocks_fmt, vbmat_entries_fmt);

        if (params.verbose > 0)    cout << "VBS matrix (Asymmetric Angle Method) created:" << endl;
        if (params.verbose > 1)    matprint(vbmat_algo);

        delete[] algo_col_part;

        //size of minimum, mazimum, average height of nonzero blocks.
        intT max_block_H = 0;
        intT min_block_H = params.A_rows;
        float avg_block_height = 0.;
        intT tot_nz_blocks = 0;
        for (intT i = 0; i < vbmat_algo.block_rows; i++)
        {
            intT b_size = vbmat_algo.row_part[i + 1] - vbmat_algo.row_part[i];
            avg_block_height += b_size * vbmat_algo.nzcount[i];
            tot_nz_blocks += vbmat_algo.nzcount[i];
            if (b_size > max_block_H) max_block_H = b_size;
            if (b_size < min_block_H) min_block_H = b_size;
        }
        avg_block_height /= tot_nz_blocks;

        //accumulate results in vectors
        avg_height_vec.push_back(avg_block_height);
        total_area_vec.push_back(vbmat_algo.nztot);
        block_rows_vec.push_back(vbmat_algo.block_rows);
        nz_blocks_vec.push_back(tot_nz_blocks);
        min_block_vec.push_back(min_block_H);
        max_block_vec.push_back(max_block_H);

        cleanVBS(vbmat_algo);
    }
    
    cleanCSR(input_cmat);


    //*******************************************
    //	 COLLECT EXPERIMENTS RESULTS
    //******************************************

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



    //OUTPUT PHASE
    if ((params.verbose == -1) or (params.verbose > 1))
    {
        cout << output_names << endl;
        cout << output_values << endl;
    }

}