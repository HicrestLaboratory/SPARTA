#include <vector>
#include "utilities.h"
#include <numeric> //std::iota

using namespace std;


vector<intT> get_permutation(const vector<intT> &grouping)
{
    auto comp = [&](int i, int j)
    {
        return grouping[i] < grouping[j];
    };

    vector<intT> v(grouping.size());
    iota(v.begin(), v.end(), 0);
    sort (v.begin(), v.end(), comp);

    return v;
}

vector<intT> get_partition(const vector<intT> &grouping)
{
    vector<intT> reordered_groups(grouping);
    vector<intT> partition;

    sort (reordered_groups.begin(), reordered_groups.end());
    
    intT current_group = -1;
    intT current_size = 0;
    intT i = 0;
    while (i < reordered_groups.size())
    {
        if (current_group != reordered_groups[i])
        {
            current_group = reordered_groups[i];
            partition.push_back(i);
        }
        i++;
    }
    partition.push_back(reordered_groups.size());
    return partition;
}
vector<intT> merge_rows(vector<intT> A, intT*B, intT size_B)
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