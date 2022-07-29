#import "blocking.h"
#import "matrices.h"
#import "similarity.h"
#import <iostream>

using namespace std;

void IterativeBlockingJaccard(const CSR& cmat, intT* grouping)
{
    float tau = 0.5; //TODO GRAB TAU FROM ENV OR OBJECT
    for (intT i = 0; i < cmat.rows; i++) grouping[i] = -1;

    for (intT i = 0; i < cmat.rows; i++)
    {
        if (grouping[i] == -1)
        {
            grouping[i] = i;
            for (intT j = i + 1; j < cmat.rows; j++)
            {
                if (grouping[j] == -1)
                {
                    float dist = JaccardDistance(cmat.ja[i], cmat.nzcount[i], cmat.ja[j], cmat.nzcount[j]);
                    if (dist < tau) grouping[j] = i;
                }
            }
        }
    }
}

void IterativeBlockingGeneral(const CSR& cmat, intT* grouping, float tau, distFunc distanceFunction)
{
    for (intT i = 0; i < cmat.rows; i++) grouping[i] = -1;

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
}

void IterativeBlockingPattern(const CSR& cmat, intT* grouping, float tau, distFuncGroup distanceFunction)
{
    for (intT i = 0; i < cmat.rows; i++) grouping[i] = -1;

    for (intT i = 0; i < cmat.rows; i++)
    {
        if (grouping[i] == -1)
        {
            vector<intT> pattern = new vector<intT>();
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