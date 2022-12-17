#include <vector>
#include "utilities.h"

using namespace std;


vector<intT> merge_rows(vector<intT> A, intT*B, intT size_B)
{
    intT size_A = A.size();
    //A,B sparse rows (compressed indices format)
    intT i = 0;
    intT j = 0;
    vector<intT> result_vec;
    while (i < size_A && j < size_B)
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