#include "cuda_utilities.h"
#include "matrices.h"
#include "blocking.h"
#include "utilities.h"
#include "definitions.h"

#include <string.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <random>
#include <algorithm> //fill, equal


using namespace std;


int main(int argc, char* argv[])
{
    //import CSR
    CLineReader cli(argc, argv);
    if (cli.verbose_ > 0) cli.print();
    CSR cmat(cli); //sparse operand

    intT A_rows = cmat_A.rows;
    intT A_cols = cmat_A.cols;
    intT B_rows = A_cols;
    intT B_cols = cli.B_cols_;
    intT C_cols = B_cols;
    intT C_rows = A_rows;

    //dense operand randomly initialized
    random_device rd;
    mt19937 e2(rd());
    uniform_real_distribution<> dist(0, 1);

    DataT* mat_B = new DataT[B_rows * B_cols];
    for (int n = 0; n < B_rows*B_cols; n++) 
    {
        mat_B[n] = dist(e2);
    }
    //______________________________________

    DataT* mat_C = new DataT[C_rows * C_cols];

    MultiplicationAlgo mult_algo = static_cast<MultiplicationAlgo>(cli.multiplication_algo_);
    //ensure compatibility of blocking with selected multiplication 
    if (mult_algo == cusparse_bellpack)
    {
        if !(cli.force_fixed_size || cli.blocking_algo_ == fixed_size) {
            cli.force_fixed_size = true;
            cout << "WARNING: forcing fixed size. -F set to 1" << endl;
        }

        if (cli.row_block_size_ != cli.col_block_size_) {
            cli.row_block_size_ = cli.col_block_size_;
            cout << "WARNING: forcing square blocks. -b set equal to -B" << endl;
        }
    }
    //_________________________________________________________________


    BlockingEngine bEngine(cli);


    //keep track of time
    float dt;
    vector<double> algo_times;
    float mean_time;
    float std_time;

    //multiply according to multiplication algo selection
    switch (mult_algo)
    {
    case NO_MULT:
        break;

    case cublas_gemm:
        //convert to dense. 
        //Run dense multiplication
        cout << "GEMM NOT IMPLEMENTED YET" << endl;
        break;
    
    case cublas_vbr:
        bEngine.GetGrouping(cmat);
        VBR vbmat_cublas;
        vbmat_cublas.fill_from_CSR_inplace(cmat, bEngine.grouping_result, cli.col_block_size_);
        algo_times.clear();
        for (int i = -cli.warmup_; i < cli.exp_repetitions_; i++)
        {
            fill(mat_C, mat_C + C_cols*C_rows, 0);
            cublas_blockmat_multiplyAB(vbmat_cublas, mat_B, B_cols, mat_C, dt, cli.n_streams_);
            if (i >= 0) algo_times.push_back(dt); //only saves non-warmup runs
        }
        bEngine.multiplication_timer_avg = avg(algo_times);
        bEngine.multiplication_timer_std = 0; //TODO add function to calculate error in utilities.
        break;

    case cusparse_spmm:
        //TODO CONVERT AND RUN MULTIPLICATION

        //bEngine.multiplication_timer_avg = avg(algo_times);
        //bEngine.multiplication_timer_std = 0; //TODO add function to calculate error in utilities.
        break;

    case cusparse_bellpack:
        bEngine.GetGrouping(cmat);
        VBR vbmat_bellpack;
        vbmat_bellpack.fill_from_CSR_inplace(cmat, bEngine.grouping_result, cli.col_block_size_);
        //TODO CONVERT AND RUN MULTIPLICATIONS
        
        //bEngine.multiplication_timer_avg = avg(algo_times);
        //bEngine.multiplication_timer_std = 0; //TODO add function to calculate error in utilities.
        break;
    }

    ofstream outfile;
    bool save_grouping = false;
    outfile.open(cli.outfile_);

    ofstream outfile_grouping;
    if (save_grouping)
    {
        outfile_grouping.open(cli.outfile_ + ".g");
    }
    save_blocking_data(outfile, cli, bEngine, cmat, save_grouping, outfile_grouping);


}
