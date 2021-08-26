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
}


template <class myType>
int output_couple(std::string& names, std::string& values, std::string name, myType value)
{
    //append name to names and value to values; add spaces
    //used to produce output in CSV-like form;

    names += name + " ";
    values += std::to_string(value) + " ";
}

int output_couple(std::string& names, std::string& values, std::string name, std::string value);


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
    std::string reorder_algo = "saad";

    int algo_block_size = 4;
    
    std::string input_source;
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

int output_couple_parameters(input_parameters&, std::string& names, std::string& values);

int get_input_params(int argc, char* argv[], input_parameters& params);

int get_input_CSR(CSR& cmat_A, input_parameters& params);

int scramble_input(CSR& cmat, input_parameters& params);



