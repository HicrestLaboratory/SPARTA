#include "matrices.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include "utilities.h"

using namespace std;

int main(int argc, char* argv[])
{
    CLineReader cli(argc, argv);
    if (cli.verbose_ > 0) cli.print();
    CSR cmat(cli);
    if (cli.verbose_ > 0) cmat.print(cli.verbose_);

    /*
    cout << "converting to vbr" << endl;
    VBR vbmat;
    vector<intT> partition{0,2,5,6,9,cmat.rows};
    vbmat.fill_from_CSR(cmat, partition, 3);
    vbmat.print();
    */


    intT B_rows = cmat.cols;
    intT B_cols = 5;
    DataT* mat_B = new DataT[B_cols*B_rows];
    fill_n(mat_B,B_cols*B_rows,1);

    cout << "B:" << endl;
    cout << B_cols << endl;
    cout << B_rows << endl;
    print_mat(mat_B,B_rows,B_cols, B_rows);


    DataT* mat_C = new DataT[B_cols*cmat.rows]{0};
    cmat.multiply(mat_B,B_cols,mat_C);


    print_mat(mat_C,cmat.rows,B_cols, cmat.rows);

    cout << "TEST COMPLETED" << endl;

}