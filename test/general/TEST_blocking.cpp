#include "matrices.h"
#include "blocking.h"
#include "similarity.h"
#include <fstream>
#include <iostream>

using namespace std;

int main()
{
    CSR cmat;
    ifstream fin;
    fin.open("data/TEST_matrix_weighted.txt");
    cmat.read_from_edgelist(fin, " ", false);
    cmat.print();
    intT* grouping = new intT[cmat.rows];


    distFunc distanceFunction = &JaccardDistance;
    float tau = 0.7;
    IterativeBlockingGeneral(cmat, grouping, tau, distanceFunction);

    cout << "MATRIX BLOCKED WITH JACCARD: BLOCKING =" << endl;
    for (intT i = 0; i < cmat.rows; i++)
    {
        cout << grouping[i] << endl;
    }

    intT block_size = 3;
    distFunc distanceFunction = [=](intT* A, intT a, intT*B,intT b){return HammingDistanceQuotient(A,a,B,b,block_size);};
    float tau = 0.7;
    IterativeBlockingGeneral(cmat, grouping, tau, distanceFunction);

    cout << "MATRIX BLOCKED WITH JACCARD: BLOCKING =" << endl;
    for (intT i = 0; i < cmat.rows; i++)
    {
        cout << grouping[i] << endl;
    }


    cmat.clean();
}