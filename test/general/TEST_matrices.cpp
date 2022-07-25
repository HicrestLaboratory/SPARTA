#include "matrices.h"
#include <fstream>
#include <iostream>

using namespace std;

int main()
{
    CSR cmat = new CSR();
    ifstream fin;
    fin.open("data/TEST_matrix_weighted.txt");
    cmat.read_from_edgelist(fin, " ", false);
    cmat.print();
}