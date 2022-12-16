#include <vector>
#include "utilities.h"

using namespace std;


vector<intT> merge_rows(intT* A, intT sizeA, intT*B, intT sizeB)
{
    //A,B sparse rows (compressed indices format)
    intT i = 0;
    intT j = 0;
    vector<intT> result_vec;
    while (i < sizeA && j < sizeB)
    {
        if (A[i] <= A[j])
        {
            i++;
            result_vec.push_back(A[i]);
        }   

        if (A[i] > A[j])
        {
            j++;
            result_vec.push_back(A[j]);
        }
    }
    return result_vec;
}