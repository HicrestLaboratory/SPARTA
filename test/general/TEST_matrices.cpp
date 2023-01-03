#include "matrices.h"
#include <fstream>
#include <iostream>

using namespace std;

int main()
{
    ifstream fin;
    fin.open("data/TEST_matrix_weighted.txt");
    CSR cmat(fin, " ", false);
    cmat.print(1);

    cout << "converting to vbr" << endl;
    VBR vbmat;
    vector<intT> partition{0,2,5,6,9};
    vbmat.fill_from_CSR(cmat, partition, 3);
    vbmat.print();
}