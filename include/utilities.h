#include <vector>
#include <algorithm> //std::copy
#include <iostream>
#include <fstream>
#include <numeric> //std::accumulate
#include "input.h"
#include "matrices.h"
#include "blocking.h"

std::vector<intT> merge_rows(std::vector<intT> A, intT*B, intT sizeB);

std::vector<intT> get_partition(const std::vector<intT> &grouping);

std::vector<intT> get_permutation(const std::vector<intT> &grouping);

std::vector<intT> get_fixed_size_grouping(const std::vector<intT> &grouping, intT row_block_size);


bool check_structured_sparsity(std::vector<intT>& structured_sparsity_pattern, std::vector<intT>& structured_sparsity_column_counter, intT* row, intT row_len, int m);
void update_structured_sparsity(std::vector<intT>& structured_sparsity_pattern, std::vector<intT>& structured_sparsity_column_counter, intT* row, intT row_len);


template<typename T>
double avg(std::vector<T> const& v) {
    if (v.empty()) 
    {
        return 0;
    }

    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}


void save_blocking_data(std::ostream &outfile, CLineReader &cLine, BlockingEngine &bEngine, CSR &cmat, bool save_blocking = false, std::ostream &blocking_outfile = std::cout);

template <class MyType>
void print_vec(std::vector<MyType> vec, std::ostream& stream = std::cout, std::string separator = " ")
{
    for (intT i = 0; i < vec.size(); i++)
    {
        stream << vec[i] << separator;
    }
    stream << std::endl;
}


template<typename T>
void print_mat(T* mat, intT rows, intT cols, intT main_dim, bool rowwise = false, std::ostream& stream = std::cout) 
{
    //print mat in rowwise or columnwise format
    for (intT i = 0; i < rows; i++)
    {
        for (intT j = 0; j < cols; j++)
        {
            T val = rowwise? mat[j + cols*i] : mat[i + rows*j];
            stream << val << " ";
        }
        stream << std::endl;
    }
    stream << std::endl;
}

//permutes an array of n elements (original) according to a permutation (perm);
template <class myType>
void permute(myType* arr, std::vector<intT> &perm) {
    //TODO efficient version with cycles
    myType* temp_arr = new myType[perm.size()];
    std::copy(arr, arr + perm.size(), temp_arr);   //auto copying

    for (intT i = 0 ; i < perm.size(); i++)
    {
        arr[i] = temp_arr[perm[i]];
    }

    delete[] temp_arr;
}