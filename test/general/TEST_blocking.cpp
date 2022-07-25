#include "matrices.h"
#include "blocking.h"
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
    saad_reordering(cmat, grouping);

    cout << "MATRIX BLOCKED: BLOCKING =" << endl;
    for (intT i = 0; i < cmat.rows; i++)
    {
        cout << grouping[i] << endl;
    }

    cmat.clean();
}