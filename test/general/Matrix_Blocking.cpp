#include "matrices.h"
#include "blocking.h"
#include "utilities.h"

#include <fstream>
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
    CLineReader cli(argc, argv);
    if (cli.verbose_ > 1) cli.print();
    CSR cmat(cli);
    if (cli.verbose_ > 1) cmat.print(cli.verbose_);
    BlockingEngine bEngine(cli);

    //evaluate the grouping
    vector<intT> grouping = bEngine.GetGrouping(cmat); 
    bEngine.print();

    ofstream outfile;
    bool save_grouping = false;
    outfile.open(cli.outfile_);


    ofstream outfile_grouping;
    if (save_grouping)
    {
        outfile_grouping.open(cli.outfile_ + ".g");
    }
    save_blocking_data(outfile, cli, bEngine, cmat, save_grouping, outfile_grouping);

    cout << "TEST COMPLETED" << endl;
}