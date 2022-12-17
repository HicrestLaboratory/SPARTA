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

    distFunc distanceFunction = &JaccardDistance;
    float tau = 0.7;
    vector<intT> grouping_1 = IterativeBlockingGeneral(cmat, tau, distanceFunction);

    cout << "MATRIX BLOCKED WITH (BASIC) JACCARD: " << endl;
    cmat.reorder(grouping_1);
    cmat.print();


    intT block_size = 3;
    tau = 0.7;
    distFuncGroup distanceFunctionGroup = JaccardDistanceGroup;
    vector<intT> grouping_2 = IterativeBlockingPattern(cmat, tau, distanceFunctionGroup, block_size);
    cmat.reorder(grouping_2);

    cmat.print();
    cmat.clean();
}