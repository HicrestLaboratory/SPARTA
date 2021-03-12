#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <unistd.h>
#include <math.h>
#include <typeinfo>
#include <iterator>
#include <algorithm>

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


int main(int argc, char* argv[]) {

    opterr = 0;

    int verbose = 3;

    intT A_rows = 12;            //rows in the square input matrix;
    intT A_cols = 8;
    intT mat_A_entries_fmt = 1;          //format inside blocks. 0: row_major; 1: col_major
    intT mat_A_block_fmt = 1;            //format of blocks. 0: compressed sparse rows. 1: compressed sparse columns. 
    intT row_block_size = 4;             //block size for variable block matrix.
    intT col_block_size = 4;
    float density = 0.5;        //density of the input matrix;
    float block_density = 0.5;  //density inside the blocks;

    intT B_cols = 5;             //number of columns in the output matrix;
    float B_density = 1.;       //density of the multiplication matrix

    int seed = 123;




    //terminal options loop
    opterr = 0;
    char c;
    while ((c = getopt(argc, argv, "b:q:Q:m:n:p:P:k:S:f:F:")) != -1)
        switch (c)
        {
        
        case 'b': //block density
            //has only effect for example 4
            block_density = stof(optarg);
            if (block_density < 0 or block_density > 1) {
                fprintf(stderr, "Option -b tried to set block density outside of [0,1]");
                return 1;
            }
            break;

        case 'm': //input matrix rows
            A_rows = stoi(optarg);
            break;

        case 'n': //output matrix cols
            B_cols = stoi(optarg);
            break;

        case 'k': //input matrix columns
            A_cols = stoi(optarg);
            break;

        case 'f': //entries format (0 or 1)
            mat_A_entries_fmt = stoi(optarg);
            break;

        case 'F': //block format (0 or 1)
            mat_A_block_fmt = stoi(optarg);
            break;

        case 'p': //size of row blocks
            row_block_size = stoi(optarg);
            break;

        case 'P': //size of col blocks
            col_block_size = stoi(optarg);
            break;

        case 'q': //density of input matrix
            density = stof(optarg);
            if (density < 0 or density > 1) {
                fprintf(stderr, "Option -k tried to set density outside of [0,1]");
                return 1;
            }
            break;

        case 'Q': //density of input matrix
            block_density = stof(optarg);
            if (block_density < 0 or block_density > 1) {
                fprintf(stderr, "Option -k tried to set block_density outside of [0,1]");
                return 1;
            }
            break;

        case 'S': //random seed
            seed = stoi(optarg);
            break;

        case 'v': //verbose
            verbose = stoi(optarg);
            break;

        case '?':
            fprintf(stderr, "Option -%c does not exists, or requires an argument.\n", optopt);
            return 1;
        default:
            abort();
        }

    //TODO -h HELP

    srand(seed);

    //GENERATE A RANDOM BLOCK SPARSE MATRIX

    //A_rows and density have been previously set by options.

    VBS vbmat;

    random_sparse_blocks_mat(vbmat, A_rows, A_cols, mat_A_block_fmt, mat_A_entries_fmt, row_block_size, col_block_size, block_density, density);
       
    if (verbose > 0)
    {
        cout << "CREATED A RANDOM BLOCK MATRIX:"
            << " Rows = " << A_rows << "\n"
            << " Columns: " << A_cols << "\n"
            << " Row Block size: " << row_block_size << "\n"
            << " Col Block size: " << col_block_size << "\n"
            << " Density OF blocks: " << block_density << "\n"
            << " Density IN blocks: " << density << "\n"
            << endl;
    }

    if (verbose > 1) //print it
    {
        matprint(vbmat);
    }

    //___________________________________________
    //*******************************************
    //		END OF INPUT
    //*******************************************
   
    //output format
    cout << fixed;

    //creating a random matrix B
    intT B_rows = A_cols;
    int mat_B_fmt = 1;
    DataT* mat_B = new DataT[B_rows * B_cols]{ 0 };
    random_mat(mat_B, B_rows, B_cols, B_density); // creates a random DataT matrix filled with 1s at a fixed density

    if (verbose > 0)        std::cout << "Random matrix B created:" << std::endl;
    if (verbose > 1)        matprint(mat_B, B_rows, B_cols, B_rows, mat_B_fmt);





    //Cleaning phase
    delete[] mat_B;
    cleanVBS(vbmat);

}
 
 