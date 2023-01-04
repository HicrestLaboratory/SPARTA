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

template<typename T>
double avg(std::vector<T> const& v) {
    if (v.empty()) 
    {
        return 0;
    }

    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}


void save_blocking_data(std::ostream &outfile, CLineReader &cLine, BlockingEngine &bEngine, CSR &cmat, bool save_blocking = false);

template <class MyType>
void print_vec(std::vector<MyType> vec, std::ostream& stream = std::cout, std::string separator = " ")
{
    for (intT i = 0; i < vec.size(); i++)
    {
        stream << vec[i] << separator;
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