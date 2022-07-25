#import "blocking.h"
#import "matrices.h"
#import "similarity.h"

void saad_reordering(const CSR& cmat, intT* grouping)
{
    float tau = 0.5; //TODO GRAB TAU FROM ENV OR OBJECT
    grouping = new intT[cmat.rows];
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
                    j = JaccardDistance(cmat.ja[i], cmat.nzcount[i], cmat.ja[j], cmta.nzcount[j]);
                    if (j < tau) grouping[j] = j;
                }
            }
        }
    }
}
