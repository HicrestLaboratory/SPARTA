#include "matrices.h"
#include "blocking.h"
#include <fstream>
#include <iostream>

using namespace std;

int main()
{

    //Handle input
    
    //Determine Blocking parameters from input
    BlockingEngine BEngine;
    BEngine.tau = 0.1;
    BEngine.block_size = 3;
    BEngine.use_groups = false;
    BEngine.use_pattern = true;
    BEngine.SetComparator(1);

    
    //import matrix
    ifstream fin;
    fin.open("data/TEST_matrix_weighted.txt");
    CSR cmat(fin, " ", false);


    //determine grouping
    cout << "evaluating reordering" <<endl;
    vector<intT> grouping = BEngine.ObtainPartition(cmat);

    //calculate statistics
    vector<intT> nz_block_count = cmat.get_VBR_nzcount(grouping, BEngine.block_size);
    vector<intT> partition = get_partition(grouping);
    intT total_nonzeros = calculate_nonzeros_VBR(partition, nz_block_count, BEngine.block_size);

    //save results on file

}