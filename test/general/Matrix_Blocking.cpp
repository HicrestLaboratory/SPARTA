#include "matrices.h"
#include "blocking.h"
#include "utilities.h"

#include <fstream>
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
    CLineReader cli(argc, argv);
    if (cli.verbose_ > 0) cli.print();
    CSR cmat(cli);
    BlockingEngine bEngine(cli);

    //evaluate the grouping
    vector<intT> grouping = bEngine.GetGrouping(cmat); 
    bEngine.print();

    ofstream outfile;
    bool save_grouping = true;
    outfile.open(cli.outfile_);
    save_blocking_data(outfile, cli, bEngine, cmat, save_grouping);

    cout << "TEST COMPLETED" << endl;
}