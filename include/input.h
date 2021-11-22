#pragma once
#include <string>
#include <vector>
#include "comp_mats.h"

typedef std::vector<float> vec_d;
typedef std::vector<std::string> vec_str;

template<class T>
float mean(std::vector<T> v)
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

template<class T>
float std_dev(std::vector<T> v)
{
    //std of a vector
    float m = mean(v);
    float s = 0.;
    for (auto t : v)
    {
        s += (t - m) * (t - m);
    }
    s = sqrt(s / v.size());
    return s;
}


template <class myType>
int output_couple(std::string& names, std::string& values, std::string name, myType value)
{
    //append name to names and value to values; add spaces
    //used to produce output in CSV-like form;

    names += name + " ";
    values += std::to_string(value) + " ";

    return 0;
}

int output_couple(std::string& names, std::string& values, std::string name, std::string value);

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
        intT min_block_H = INT_MAX;
        float avg_block_height = 0.;
        intT tot_nz_blocks = 0;
        for (intT i = 0; i < vbmat.block_rows; i++)
        {
            intT b_size = vbmat.row_part[i + 1] - vbmat.row_part[i];
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

struct input_parameters
{
    int verbose = 3;
    int input_type = 4;

    intT A_rows = 12;            //rows in the square input matrix;
    intT A_cols = 8;
    int cmat_A_fmt = 0;          //cuda needs column-major matrices
    int A_nnz = -1;              //will store the number of nonzeros

    
    intT B_cols = 5;             //number of columns in the output matrix;
    float B_density = 1.;       //density of the multiplication matrix
    
    int block_size = 4;         //block size for variable block matrix. Rows and columns must be evenly divisible by this;
    float density = 0.5;        //density of the input matrix;
    float block_density = 0.5;  //density inside the blocks;

    std::string exp_name = "default";
    std::string reorder_algo = "saad_blocks";
    std::string similarity_func = "scalar";
    int hierarchic_merge= 1;         //Activate hierchical merging?
    float merge_limit = -1;            //the merge limit. If -1, use the theoretical limit; if 0, deactivate;


    int algo_block_size = 4;
    bool save_reordering = false; //save the grouping and blocking info obtained from reordering?

    std::string input_source = "NO_INPUT";
    int scramble = 0;           //scramble the input matrix?
    int n_streams = 16;         //number of streams for the custom block_based multiplication

    float eps = 0.5;            //this value sets how different two rows in the same block can be.
                                //eps = 1 means only rows with equal structure are merged into a block
                                //eps = 0 means all rows are merged into a single block
    int seed = 123;
    float precision = 0.0001;   //precision for float equality check

    int warmup = 0;             //number of warmup experiments
    int experiment_reps = 5;    //number of non-warmup repetitions
    int algo = -1;              //algorithm choice (-1: all)
    int check_correct = 0;      //verify correctness?

};

int output_couple_parameters(input_parameters& params, std::string& names, std::string& values);

int save_reordering(std::string& output_file, intT* hash_groups, input_parameters& params);

int get_input_params(int argc, char* argv[], input_parameters& params);

int get_input_CSR(CSR& cmat_A, input_parameters& params);

int scramble_input(CSR& cmat, input_parameters& params);



