#include <vector>
#include "utilities.h"

using namespace std;


std::vector<intT> merge_rows(std::vector<intT> A, intT*B, intT size_B)
{
    //A,B sparse rows (compressed indices format)
    intT i = 0;
    intT j = 0;
    vector<intT> result;
    intT size_A = A.size();

    while (i < size_A && j < size_B)
    {
        if (A[i] <= B[j])
        {
            i++;
            result.push_back(A[i]);
        }   

        if (A[i] > B[j])
        {
            j++;
            result.push_back(B[j]);
        }
    }
    return result;
}