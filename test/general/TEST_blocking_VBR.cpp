#include "matrices.h"
#include "blocking.h"
#include "utilities.h"

#include <fstream>
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
    CLineReader cli(argc, argv);

    CSR cmat(cli);
    cmat.print();

    BlockingEngine BEngine(cli);

    //evaluate the grouping
    cout << "evaluating reordering" <<endl;
    vector<intT> grouping = BEngine.ObtainPartition(cmat);


    cout << "create VBR from grouping;" << endl;
    //create a VBR matrix from grouping (without reordering the original csr)
    VBR vbmat;
    vbmat.fill_from_CSR_inplace(cmat, grouping, 3);
    vbmat.print();
    
    vector<intT> nzcount_VBR = cmat.get_VBR_nzcount(grouping,3);
    cout << "nzcount_VBR: "; 
    print_vec(nzcount_VBR);

    cout << "create VBR from row_partition;" << endl;
    //create a VBR matrix from the row_partition (reordering the original csr)
    VBR vbmat2; 
    cmat.reorder(grouping);
    cmat.print();
    vbmat2.fill_from_CSR(cmat, get_partition(grouping), 3);
    vbmat2.print();
}