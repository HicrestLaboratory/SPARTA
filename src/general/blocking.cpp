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
                    cout << i << " " << j << " " << dist << endl;
                }
            }
        }
    }
}