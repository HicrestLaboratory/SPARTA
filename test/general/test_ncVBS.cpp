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
#include "nvbs_utilities.h"
#include "comp_mats.h"

int main(int argc, char* argv[]) {

    //input
    {
        int verbose = 3;

        //matrix features
        intT mat_rows = 12;
        intT mat_cols = 8;
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
        intT algo_block_size = 3;
        int experiment_reps = 1;
        int scramble = 3;
        string exp_name = "default";

        //terminal options loop
        int opterr = 0;
        char c;
        while ((c = getopt(argc, argv, "b:e:f:i:m:n:p:P:q:r:S:s:v:")) != -1)
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
                        //                          1: SCRAMBLE the rows 
                        //                          2: SCRAMBLE the columns
                        //                          3: SCRAMBLE both rows and columsn
                scramble = stoi(optarg);
                break;

            case 'v': //verbose
                verbose = stoi(optarg);
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
    }


    ncVBS vbmat;
    int mat_rows = 10;
    int mat_cols = 20;
    int block_size = 5;
    float mat_density = 0.1f;
    float row_density = 0.2f;

    random_ncVBS(vbmat, mat_rows, mat_cols, block_size, mat_density, row_density)




}