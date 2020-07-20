#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>


#include <vector>
#include <string>


#include "sparse_utilities.h"
#include "comp_mats.h"

using namespace std;
typedef std::vector<float> vec_d;
typedef std::vector<string> vec_str;
typedef std::vector<int> vec_i;

float mean(vec_d v)
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

float std_dev(vec_d v)
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

int main(int argc, char* argv[]) {


    int verbose = 3;
    
    //matrix features
    int mat_rows = 12;
    int mat_cols = 8;
    int mat_fmt = 0;

    string input_source;
    int seed = 123;
    int input_block_size = 4;
    int input_block_density = 0.5;
    int input_entries_density = 0.5;
    int* A_row_part;
    int* A_col_part;

    //algorithm and experiments parameters
    float eps;
    int generate_new_random = 0;

    int algo_block_size = 6;
    int experiment_rep = 5;
    int scramble = 1; //scramble the input matrix?

    //terminal options loop
    int opterr = 0;
    char c;
    while ((c = getopt(argc, argv, "b:e:f:i:m:n:p:P:q:r:S:s")) != -1)
        switch (c)
        {
        case 'i':// select input example
            input_type = stoi(optarg);
            //  1: Random CSR
            //  2: SNAP Edgelist
            //  3: MTX Format
            //  4: Random Variable Block matrix
            break;

        case 'b': //density of blocks (% of nonzero blocks)
            //has only effect for example 4
            block_density = stof(optarg);
            if (block_density < 0 or block_density > 1) {
                fprintf(stderr, "Option -b tried to set block density outside of [0,1]");
                return 1;
            }
            break;

        case 'q': //density of entries. if i = 4, this is the density INSIDE each block.
            //has only effect for example 1 and 4
            density = stof(optarg);
            if (density < 0 or density > 1) {
                fprintf(stderr, "Option -k tried to set density outside of [0,1]");
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
            mat_rows = stoi(optarg);
            break;

        case 'n': //input matrix cols
            //has only effect for example 1 and 4
            mat_cols = stoi(optarg);
            break;

        case 'p': //size of input blocks
            //ony used if i = 4, random VBS
            input_block_size = stoi(optarg);
            break;

        case 'P': //size of the blocks used by the reordering algorithm
            algo_block_size = stoi(optarg);
            break;

        case 'r': //number of experiment repetitions
            experiment_reps = stoi(optarg);
            break;

        case 's': //scramble the input matrix?  0: NO SCRAMBLE
                    //                          1: SCRAMBLE the rows randomly 
            scramble = stoi(optarg);
            break;

        case 'S': //random seed; use -1 for random seeding.
            seed = stoi(optarg);
            if (seed == -1) seed = time(NULL);
            srand(seed);
            break;

        case 'v': //verbose
            verbose = stoi(optarg);
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

    for (int current_repetition = 0; current_repetition < experiment_rep; current_repetition++)
    {
        //INPUT EXAMPLE 1: RANDOM CSR
        //create a random sparse matrix
        if (input_type == 1) {
            if ((current_repetition == 0) or (generate_new_random == 1))
            {
                DataT* rand_mat = new DataT[mat_cols * mat_rows];
                random_mat(rand_mat, mat_rows, mat_cols, input_entries_density); //generate random mat //todo: generate directly the CSR

                convert_to_CSR(rand_mat, mat_rows, mat_cols, mat_fmt, input_cmat, input_cmat_fmt);
                delete[] rand_mat;

                if (verbose > 0) cout << "CREATED A RANDOM CSR with density = " << density << endl;

            }
        }
        //______________________________________


        //TEST
        //INPUT EXAMPLE 2: read graph in edgelist format into CSR
        if (input_type == 2) {
            if (current_repetition == 0)
            {
                if (input_source.empty()) input_source = "testgraph1.txt";

                string delimiter = " ";
                GraphMap snap_graph;
                read_snap_format(snap_graph, input_source, delimiter);         //Read into a GraphMap matrix from a .txt edgelist (snap format)
                MakeProper(snap_graph);
                convert_to_CSR(snap_graph, input_cmat, input_cmat_fmt);

                if (verbose > 0) cout << "IMPORTED A CSR FROM A SNAP EDGELIST" << endl;
            }

        }
        //______________________________________


        //INPUT EXAMPLE 3: read from MTX format
        if (input_type == 3)
        {
            if (current_repetition == 0)
            {
                //read from mtx
                if (input_source.empty()) input_source = "testmat.mtx";
                read_mtx_format(input_cmat, input_source, input_cmat_fmt); //read into CSR

                if (verbose > 0)            cout << "IMPORTED A CSR FROM MTX FILE" << endl;
            }
        }
        //______________________________________


        //INPUT EXAMPLE 4: create a random matrix with block structure
        if (input_type == 4) {
            if ((current_repetition == 0) or (generate_new_random == 1))
            {
                DataT* rand_block_mat = new DataT[A_rows * A_cols];

                //TODO do not start by array but create directly the CSR?
                //TODO arbitrary partition

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
        }
        //___________________________________________


        //*******************************************
        //		END OF INPUT
        //spmat must hold a proper CSR matrix at this point
        //******************************************

        if (verbose > 0) cout << "INPUT ACQUIRED." << endl;
        if (verbose > 1) matprint(cmat_A);

        //update rows and cols count to input values
        mat_rows = input_cmat.rows;
        mat_cols = input_cmat.cols;

        //scramble the original matrix
        if (scramble)
        {

            if (verbose > 0) cout << "input matrix rows scrambled" << endl;
            int* random_permutation = new int[mat_rows];
            randperm(random_permutation, A_rows);
            permute_CSR(input_cmat, random_permutation, 0);
            delete[] random_permutation;
        }


    }
}