#include "blocking.h"
#include "matrices.h"
#include "similarity.h"
#include <vector>
#include <iostream>

using namespace std;

std::vector<intT> IterativeBlockingGeneral(const CSR& cmat, float tau, distFunc distanceFunction)
{
    vector<intT> grouping(cmat.rows, -1);

    for (intT i = 0; i < cmat.rows; i++)
    {
        if (grouping[i] == -1)
        {
            grouping[i] = i;
            for (intT j = i + 1; j < cmat.rows; j++)
            {
                if (grouping[j] == -1)
                {
                    float dist = distanceFunction(cmat.ja[i], cmat.nzcount[i], cmat.ja[j], cmat.nzcount[j]);
                    if (dist < tau) grouping[j] = i;
                }
            }
        }
    }
    return grouping;
}


/*
void IterativeBlockingPattern(const CSR& cmat, intT* grouping, float tau, distFuncGroup distanceFunction)
{
    for (intT i = 0; i < cmat.rows; i++) grouping[i] = -1;

    for (intT i = 0; i < cmat.rows; i++)
    {
        if (grouping[i] == -1)
        {
            vector<intT> pattern;
            intT current_group_size = 1;
            grouping[i] = i;
            pattern.insert(pattern.end(), &cmat.ja[i][0], &cmat.ja[i][cmat.nzcount[i]]); //copy the current row into the pattern

            for (intT j = i + 1; j < cmat.rows; j++)
            {
                if (grouping[j] == -1)
                {
                    float dist = distanceFunction(pattern, pattern.size(), current_group_size, cmat.ja[j], cmat.nzcount[j], 1);
                    if (dist < tau)
                    {
                        grouping[j] = i;
                        pattern = merge_rows(pattern, pattern.size(), cmat.ja[j], cmat.nzcount[j]); //update the pattern
                        current_group_size++;
                    }
                }
            }
        }
    }
}

*/