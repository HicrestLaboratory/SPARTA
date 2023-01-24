#include "matrices.h"
#include "blocking.h"
#include <fstream>
#include <iostream>

using namespace std;

int main()
{
    ifstream fin;
    fin.open("data/TEST_matrix_weighted.txt");
    CSR cmat(fin, " ", false);
    cmat.print(1);

    BlockingEngine BEngine;
    BEngine.tau = 0.1;
    BEngine.block_size = 3;
    BEngine.use_groups = false;
    BEngine.use_pattern = true;
    BEngine.SetComparator(1);
    cout << "evaluating reordering" <<endl;
    vector<intT> grouping = BEngine.GetGrouping(cmat);

    cmat.reorder(grouping);
    cmat.print(1);


    cout << "TEST COMPLETED" << endl;

}