#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <unistd.h> //getopt, optarg
#include <ctime>
#include <cmath>
#include <algorithm> //min_element, max_element




#include <vector>
#include <string>


#include "sparse_utilities.h"
#include "reorderings.h"

#include "comp_mats.h"

using namespace std;
typedef std::vector<float> vec_d;
typedef std::vector<string> vec_str;


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
}


template <class myType>
int output_couple(string& names, string& values, string name, myType value)
{
    //append name to names and value to values; add spaces
    //used to produce output in CSV-like form;

    names += name + " ";
    values += to_string(value) + " ";
}

int output_couple(string& names, string& values, string name, string value)
{
    //append name to names and value to values; add spaces
    //used to produce output in CSV-like form;

    names += name + " ";
    values += value + " ";
}



int main(int argc, char* argv[]) {


    int verbose = 3;
    
    //matrix features
    intT mat_rows = 20;
    intT mat_cols = 40;
    intT mat_fmt = 0;

    string input_source = "no_input_source";
    int seed = 123;
    intT input_block_size = 4;
    float input_block_density = 0.5;
    float input_entries_density = 0.5;
    intT* A_row_part;
    intT* A_col_part;

    //algorithm and experiments parameters
    float eps = 0.5;
    int generate_new_random = 0;

    int input_type = 4;
    intT algo_block_size = 5;
    int experiment_reps= 1;
    int scramble = 3;
    string exp_name = "default";
    string reorder_algo = "saad";

    //terminal options loop
    int opterr = 0;
    char c;
    while ((c = getopt(argc, argv, "b:e:f:i:m:n:p:P:q:r:R:S:s:v:")) != -1)
        switch (c)
        {
        case 'N': //name of the experiment
            exp_name = optarg;

        case 'i':// select input example
            input_type = stoi(optarg);
            //  1: Random CSR
            //  2: SNAP Edgelist
            //  3: MTX Format
            //  4: Random Variable Block matrix
            break;

        case 'b': //density of blocks (% of nonzero blocks)
            //has only effect for example 4
            input_block_density = stof(optarg);
            if (input_block_density < 0 or input_block_density > 1) {
                fprintf(stderr, "Option -b tried to set block density outside of [0,1]");
                return 1;
            }
            break;

        case 'q': //density of entries. if i = 4, this is the density INSIDE each block.
            //has only effect for example 1 and 4
            input_entries_density = stof(optarg);
            if (input_entries_density < 0 or input_entries_density > 1) {
                fprintf(stderr, "Option -k tried to set entries density outside of [0,1]");
                return 1;
            }
            break;

        case 'e': //epsilon used for matrix reordering;
            eps = stof(optarg);
            if (eps < 0. or eps > 1.) {
                fprintf(stderr, "Option -e tried to set epsilon outside of [0,1]");
                return 1;
            }
            break;

        case 'f': //select source file
            //has only effect for example 2 and 3;
            input_source = optarg;
            break;

        case 'm': //input matrix rows
            //has only effect for example 1 and 4
            try
            {
                mat_rows = stoi(optarg);
                if (mat_rows <= 0) throw 1;

            }
            catch(...)
            {
                std::cout << "m (rows) should be an integer > 0";
                return 1;
            }
            break;

        case 'n': //input matrix cols
            //has only effect for example 1 and 4
            try
            {
              mat_cols = stoi(optarg);
              if (mat_cols <= 0) throw 1;

            }
            catch (...)
            {
                std::cout << "n (cols) should be an integer > 0";
                return 1;
            }
            break;
            break;

        case 'p': //size of input blocks
            //ony used if i = 4, random VBS
            try
            {
                input_block_size = stoi(optarg);
                if (input_block_size <= 0) throw 1;
            }
            catch (...)
            {
                std::cout << "p (input block size) should be an integer > 0";
                return 1;
            }
            break;

        case 'P': //size of the blocks used by the reordering algorithm
            try 
            {
                algo_block_size = stoi(optarg);
                if (algo_block_size <= 0) throw 1;
            }
            catch (...)
            {
                std::cout << "P (algorithm block size) should be an integer > 0";
                return 1;
            }
            break;

        case 'r': //number of experiment repetitions
            try
            {
                experiment_reps = stoi(optarg);
                if (experiment_reps <= 0) throw 1;
            }
            catch (...)
            {
                std::cout << "r (number of repetitions) should be an integer > 0";
                return 1;
            }
            break;

        case 's': //scramble the input matrix?  0: NO SCRAMBLE
                    //                          1: SCRAMBLE the rows 
                    //                          2: SCRAMBLE the columns
                    //                          3: SCRAMBLE both rows and columsn
            try
            {
                scramble = stoi(optarg);
                if (scramble < 0 || scramble > 3) throw 1;
            }
            catch (...)
            {
                std::cout << "s (scramble) should be an 1,2 or 3";
                return 1;
            }
            break;
        
        case 'v': //verbose
            verbose = stoi(optarg);
            break;

        case 'R': //the reordering algorithm to be used:
            //saad
            //saad-blocks
            reorder_algo = optarg;
            break;

        case 'S': //random seed; use -1 for random seeding.
            seed = stoi(optarg);
            if (seed == -1) seed = time(NULL);
            srand(seed);
            break;

        case '?':
            fprintf(stderr, "Option -%c does not exists, or requires an argument.\n", optopt);
            return 1;
        default:
            abort();


            //TODO -h HELP

        }

    //INPUT CONVERSION TO Compressed Sparse Row (CSR)

    CSR input_cmat; //this will hold the CSR matrix
    int input_cmat_fmt = 0;
    if (input_type == 1)
    {
        //INPUT EXAMPLE 1: RANDOM CSR
        //create a random sparse matrix
        DataT* rand_mat = new DataT[mat_cols * mat_rows];
        random_mat(rand_mat, mat_rows, mat_cols, input_entries_density); //generate random mat //todo: generate directly the CSR
        convert_to_CSR(rand_mat, mat_rows, mat_cols, mat_fmt, input_cmat, input_cmat_fmt);
        delete[] rand_mat;

        if (verbose > 0) cout << "CREATED A RANDOM CSR with density = " << input_entries_density << endl;
    }
    //______________________________________

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

    //INPUT EXAMPLE 4: create a random matrix with block structure
    else if (input_type == 4)
    {
        DataT* rand_block_mat = new DataT[mat_rows * mat_cols];

        //TODO arbitrary partition


        if (mat_rows % input_block_size or mat_cols % input_block_size)
        {
            //TODO exception
            std::cout << "ERROR when creating a random-sparse-blocks matrix: \n matrix dimensions (currently"
                << mat_rows << " x " << mat_cols
                << " must be a multiple of block size (" << input_block_size << ")"
                << std::endl;
            return 1;
        }


        //TODO generate directly VBS
        random_sparse_blocks_mat(rand_block_mat, mat_rows, mat_cols, mat_fmt, input_block_size, input_block_density, input_entries_density);

        convert_to_CSR(rand_block_mat, mat_rows, mat_cols, mat_fmt, input_cmat, input_cmat_fmt);

        delete[] rand_block_mat;

        if (verbose > 0)
        {
            cout << "CREATED A RANDOM BLOCK MATRIX:"
                << " Rows = " << mat_rows << "\n"
                << " Columns: " << mat_cols << "\n"
                << " Block size: " << input_block_size << "\n"
                << " Density OF blocks: " << input_block_density << "\n"
                << " Density IN blocks: " << input_entries_density << "\n"
                << endl;
        }
    }
        //___________________________________________


    //*******************************************
    //		END OF INPUT
    //spmat must hold a proper CSR matrix at this point
    //******************************************

    if (verbose > 0) cout << "INPUT ACQUIRED." << endl;
    if (verbose > 1) matprint(input_cmat);

    //update rows and cols count to input values
    mat_rows = input_cmat.rows;
    mat_cols = input_cmat.cols;


    std::cout << "rows: " << mat_rows << ", cols: " << mat_cols << std::endl;
    //*******************************************
    //	 CREATES (several) VBS from the CSR, EXTRACT FEATURES
    //******************************************

    string output_names;
    string output_values;


    intT mat_nnz = count_nnz(input_cmat);

    output_couple(output_names, output_values, "exp_name", exp_name);
    output_couple(output_names, output_values, "input_type", input_type);
    output_couple(output_names, output_values, "input_source", input_source);
    output_couple(output_names, output_values, "reorder_algorithm", reorder_algo);
    output_couple(output_names, output_values, "rows", input_cmat.rows);
    output_couple(output_names, output_values, "cols", input_cmat.cols);
    output_couple(output_names, output_values, "total_nonzeros", mat_nnz);
    output_couple(output_names, output_values, "input_blocks_density", input_block_density);
    output_couple(output_names, output_values, "input_entries_density", input_entries_density);
    output_couple(output_names, output_values, "input_block_size", input_block_size);
    output_couple(output_names, output_values, "algo_block_size", algo_block_size);
    output_couple(output_names, output_values, "epsilon", eps);

    svi total_area_vec;
    svi block_rows_vec;
    svi nz_blocks_vec;
    svi min_block_vec;
    svi max_block_vec;
    vec_d avg_height_vec;
    int scramble_cols = (scramble == 2 or scramble == 3) ? 1 : 0;
    int scramble_rows = (scramble == 1 or scramble == 3) ? 1 : 0;

    output_couple(output_names, output_values, "scramble", scramble);
    for (int current_repetition = 0; current_repetition < experiment_reps; current_repetition++)
    {
        //scramble the original matrix
        if (scramble_rows)
        {

            if (verbose > 0) cout << "input matrix rows scrambled" << endl;
            intT* random_rows_permutation = new intT[mat_rows];
            randperm(random_rows_permutation, mat_rows);
            permute_CSR(input_cmat, random_rows_permutation, 0);
            delete[] random_rows_permutation;
        }
        if (scramble_cols)
        {

            if (verbose > 0) cout << "input matrix cols scrambled" << endl;
            intT* random_cols_permutation = new intT[mat_cols];
            randperm(random_cols_permutation, mat_cols);
            permute_CSR(input_cmat, random_cols_permutation, 1);
            delete[] random_cols_permutation;
        }
        if (verbose > 1) matprint(input_cmat);



        int vbmat_blocks_fmt = 0;
        int vbmat_entries_fmt = 0;
        intT algo_block_cols = std::ceil((float)mat_cols / algo_block_size);

        //prepare the column partition
        intT* algo_col_part = new intT[algo_block_cols + 1]; 
        partition(algo_col_part, 0, input_cmat.cols, algo_block_size);

        //run the reordering algo
        reorder_parameters params;
        params.tau = eps;
        params.block_size = algo_block_size;
        params.algo = reorder_algo;
        intT* hash_groups = new intT[input_cmat.rows];
        saad_reordering(input_cmat, params, hash_groups);

        //create the block matrix
        VBS vbmat_algo;
        group_to_VBS(input_cmat, hash_groups, algo_col_part, algo_block_cols, vbmat_algo, vbmat_blocks_fmt, vbmat_entries_fmt);

        if (verbose > 0)    cout << "VBS matrix (Asymmetric Angle Method) created:" << endl;
        if (verbose > 1)    matprint(vbmat_algo);

        delete[] algo_col_part;

        //size of minimum, mazimum, average height of nonzero blocks.
        intT max_block_H = 0;
        intT min_block_H = mat_rows;
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
    if ((verbose == -1) or (verbose > 1))
    {
        cout << output_names << endl;
        cout << output_values << endl;
    }

}