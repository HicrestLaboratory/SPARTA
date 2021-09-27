#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unistd.h> //getopt, optarg

#include "input.h"
#include "comp_mats.h"
#include "sparse_utilities.h"
#include "reorderings.h"

using namespace std;

int output_couple(string& names, string& values, string name, string value)
{
    //append name to names and value to values; add spaces
    //used to produce output in CSV-like form;

    names += name + " ";
    values += value + " ";
}

int output_couple_parameters(input_parameters& params, string& output_names, string& output_values)
{
    cout << fixed;

    output_couple(output_names, output_values, "exp_name", params.exp_name);
    output_couple(output_names, output_values, "input_type", params.input_type);
    output_couple(output_names, output_values, "input_source", params.input_source);

    output_couple(output_names, output_values, "rows", params.A_rows);
    output_couple(output_names, output_values, "cols", params.A_cols);
    output_couple(output_names, output_values, "B_cols", params.B_cols);
    output_couple(output_names, output_values, "B_density", params.B_density);
    output_couple(output_names, output_values, "total_nonzeros", params.A_nnz);
    output_couple(output_names, output_values, "input_blocks_density", params.block_density);
    output_couple(output_names, output_values, "input_entries_density", params.density);
    output_couple(output_names, output_values, "input_block_size", params.block_size);

    output_couple(output_names, output_values, "reorder_algorithm", params.reorder_algo);
    output_couple(output_names, output_values, "algo_block_size", params.algo_block_size);
    output_couple(output_names, output_values, "epsilon", params.eps);
    output_couple(output_names, output_values, "similarity_func", params.similarity_func);
    output_couple(output_names, output_values, "hierarchic_merge", params.hierarchic_merge);

    output_couple(output_names, output_values, "scramble", params.scramble);
    output_couple(output_names, output_values, "Warmup", params.warmup);
    output_couple(output_names, output_values, "Repetitions", params.experiment_reps);
    output_couple(output_names, output_values, "Algorithm", params.algo);
}

int get_input_params(int argc, char* argv[], input_parameters& params)
{
    //terminal options loop
    opterr = 0;
    char c;
    while ((c = getopt(argc, argv, "a:b:B:c:i:q:e:f:F:m:M:n:p:P:r:R:k:s:S:v:w:z")) != -1)
        switch (c)
        {
        case 'i':// select input example
            params.input_type = stoi(optarg);
            //  1: Random CSR
            //  2: SNAP Edgelist
            //  3: MTX
            //  4: Random Variable Block matrix
            break;

        case 'a': //algorithm selection:
                  /*    -1: all
                         1: gemm
                         2: VBS
                         3: VBS - no zeros
                         4: VBS - AH ******DEACTIVATED*******
                         5: cusparse
                  */
            params.algo = stoi(optarg);
            break;

        case 'b': //block density
            //has only effect for example 4
            params.block_density = stof(optarg);
            if (params.block_density < 0 or params.block_density > 1) {
                fprintf(stderr, "Option -b tried to set block density outside of [0,1]");
                return 1;
            }
            break;

        case 'B': //Save reordering and blocking info
            if (stoi(optarg) == 0) params.save_reordering = false;
            else params.save_reordering = true;
            break;

        case 'c': //check correctness;
            params.check_correct = stoi(optarg);
            break;

        case 'e': //epsilon used for matrix reordering;
            params.eps = stof(optarg);
            if (params.eps < 0. or params.eps > 1.) {
                fprintf(stderr, "Option -e tried to set epsilon outside of [0,1]");
                return 1;
            }
            break;

        case 'f': //select source file
            //has only effect for example 2 and 3;
            params.input_source = optarg;
            break;

        case 'F': //similarity function
            //ALLOWED:  -hamming (distance)
            //          -jaccard
            //          -scalar
            //          -density

            params.similarity_func = optarg;
            break;

        case 'm': //input matrix rows
            //has only effect for example 1 and 4
            params.A_rows = stoi(optarg);
            break;

        case 'M': //hierarchical merge
            params.hierarchic_merge = stoi(optarg);
            break;

        case 'n': //input matrix rows
            //has only effect for example 1 and 4
            params.B_cols = stoi(optarg);
            break;

        case 'N': //name of the experiment
            params.exp_name = optarg;

        case 'k': //input matrix rows
            //has only effect for example 1 and 4
            params.A_cols = stoi(optarg);
            break;

        case 'p': //size of blocks
            //ony used if i = 4, random VBS
            params.block_size = stoi(optarg);
            break;

        case 'P': //size of the blocks used by the reordering algorithm
            try
            {
                params.algo_block_size = stoi(optarg);
                if (params.algo_block_size <= 0) throw 1;
            }
            catch (...)
            {
                std::cout << "P (algorithm block size) should be an integer > 0";
                return 1;
            }
            break;

        case 'q': //density of input matrix
            //has only effect for example 1 and 4
            params.density = stof(optarg);
            if (params.density < 0 or params.density > 1) {
                fprintf(stderr, "Option -k tried to set density outside of [0,1]");
                return 1;
            }
            break;

        case 'r': //number of experiment repetitions
            try
            {
                params.experiment_reps = stoi(optarg);
                if (params.experiment_reps <= 0) throw 1;
            }
            catch (...)
            {
                std::cout << "r (number of repetitions) should be an integer > 0";
                return 1;
            }
            break;

        case 'R': //the reordering algorithm to be used:
                  //saad
                  //saad_blocks
            params.reorder_algo = optarg;
            break;


        case 's': //scramble the input matrix?  0: NO SCRAMBLE
                    //                          1: SCRAMBLE the rows 
                    //                          2: SCRAMBLE the columns
                    //                          3: SCRAMBLE both rows and columsn
            try
            {
                params.scramble = stoi(optarg);
                if (params.scramble < 0 || params.scramble > 3) throw 1;
            }
            catch (...)
            {
                std::cout << "s (scramble) should be an 1,2 or 3";
                return 1;
            }
            break;


        case 'S': //random seed
            params.seed = stoi(optarg);
            break;

        case 'v': //verbose
            params.verbose = stoi(optarg);
            break;

        case 'w': //warmup repetitions
            params.warmup = stoi(optarg);
            break;

        case 'z': //streams for cublas_block_multiply
            params.n_streams = stoi(optarg);
            break;


        case '?':
            fprintf(stderr, "Option -%c does not exists, or requires an argument.\n", optopt);
            return 1;
        default:
            abort();
        }

        if (params.reorder_algo == "saad") params.algo_block_size = 1;

        return 0;

}

int get_input_CSR(CSR& cmat_A, input_parameters& params)
{
    //INPUT1: RANDOM CSR
    //create a random sparse matrix
    if (params.input_type == 1) {
        DataT* rand_mat = new DataT[params.A_cols * params.A_rows];
        random_mat(rand_mat, params.A_rows, params.A_cols, params.density); //generate random mat

        convert_to_CSR(rand_mat, params.A_rows, params.A_cols, params.cmat_A_fmt, cmat_A, params.cmat_A_fmt);
        delete[] rand_mat;

        if (params.verbose > 0) cout << "CREATED A RANDOM CSR with density = " << params.density << endl;
    }
    //______________________________________


    //INPUT EXAMPLE 2: read graph in edgelist format into CSR
    else if (params.input_type == 2)
    {
        //TEST
        //INPUT EXAMPLE 2: read graph in edgelist format into CSR
        if (params.input_source.empty()) params.input_source = "testgraph1.txt";

        string delimiter = "\t";
        read_edgelist(params.input_source, cmat_A, params.cmat_A_fmt, delimiter);
        if (params.verbose > 0) cout << "IMPORTED A CSR FROM A SNAP EDGELIST" << endl;

        params.A_rows = cmat_A.rows;
        params.A_cols = cmat_A.cols;
        //______________________________________
    }


    //INPUT EXAMPLE 3: read graph in MTX edgelist format into CSR
    else if (params.input_type == 3)
    {
        //INPUT EXAMPLE 3: read graph in MTX edgelist format into CSR
        if (params.input_source.empty()) params.input_source = "testgraph1.txt";

        string delimiter = " ";
        read_edgelist(params.input_source, cmat_A, params.cmat_A_fmt, delimiter);
        if (params.verbose > 0) cout << "IMPORTED A CSR FROM A SNAP EDGELIST" << endl;

        params.A_rows = cmat_A.rows;
        params.A_cols = cmat_A.cols;
        //______________________________________
    }


    //INPUT EXAMPLE 4: create a random matrix with block structure
    else if (params.input_type == 4) {


        VBS vbmat_input;
        random_sparse_blocks_mat(vbmat_input, params.A_rows, params.A_cols, 1, 1, params.block_size, params.block_size, params.block_density, params.density);

        convert_to_CSR(vbmat_input, cmat_A, 0);

        if (params.verbose > 0)
        {
            cout << "CREATED A RANDOM BLOCK MATRIX:"
                << " Rows = " << params.A_rows << "\n"
                << " Columns: " << params.A_cols << "\n"
                << " Block size: " << params.block_size << "\n"
                << " Density OF blocks: " << params.block_density << "\n"
                << " Density IN blocks: " << params.density << "\n"
                << endl;
        }

        cleanVBS(vbmat_input);

    }

    //___________________________________________
    //*******************************************
    //		END OF INPUT
    //cmat must hold a proper CSR matrix at this point
    //******************************************

    params.A_nnz = count_nnz(cmat_A);
    params.A_rows = cmat_A.rows;
    params.A_cols = cmat_A.cols;

    return 0;
}

int scramble_input(CSR& cmat, input_parameters& params)
{

    int scramble = params.scramble;
    int scramble_cols = (scramble == 2 or scramble == 3) ? 1 : 0;
    int scramble_rows = (scramble == 1 or scramble == 3) ? 1 : 0;
    
    if (scramble_rows)
    {
        if (params.verbose > 0) cout << "input matrix rows scrambled" << endl;
        intT* random_rows_permutation = new intT[params.A_rows];
        randperm(random_rows_permutation, params.A_rows);
        permute_CSR(cmat, random_rows_permutation, 0);
        delete[] random_rows_permutation;
    }
    if (scramble_cols)
    {

        if (params.verbose > 0) cout << "input matrix cols scrambled" << endl;
        intT* random_cols_permutation = new intT[params.A_cols];
        randperm(random_cols_permutation, params.A_cols);
        permute_CSR(cmat, random_cols_permutation, 1);
        delete[] random_cols_permutation;
    }
}

int save_reordering(std::string& output_file, intT* hash_groups, input_parameters& params)
{
    ofstream myfile;
    myfile.open(output_file);
    myfile << "Saved reordering\n";
    myfyle << "input_source: " << params.input_source << "\n";
    myfyle << "rows: " << params.A_rows << "\n";
    myfyle << "cols: " << params.A_cols << "\n";
    myfyle << "algo_block_size: " << params.algo_block_size << "\n";
    myfyle << "reorder_algo: " << params.reorder_algo << "\n";
    myfyle << "similarity_func: " << params.similarity_func << "\n";
    myfyle << "eps: " << params.eps << "\n";
    myfile << "grouping:";
    for (intT i = 0; i < params.A_rows; i++)
    {
        myfile << " " << hash_groups[i];
    }
    myfile << "\n";
    myfile.close();
}
