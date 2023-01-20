#include "cuda_utilities.h"
#include "matrices.h"
#include "blocking.h"
#include "utilities.h"

using namespace std;



int main(int argc, char* argv[])
{
    //keep track of time
    float dt;
    vector<double> algo_times;
    float mean_time;
    float std_time;


    //import CSR
    CLineReader cli(argc, argv);
    if (cli.verbose_ > 0) cli.print();
    CSR cmat_A(cli); //sparse operand
    BlockingEngine bEngine(cli);


    //partition for the row of cmat_A
    vector<intT> partition;
    intT part_block_size = 5;
    for (intT i = 0; i < cmat_A.rows(); i += part_block_size) partition.push_back(i);
    partition.push_back(cmat_A.rows());


    //create vbamt from equally spaced row partition of cmat_A;
    VBR vbmat; 
    if (cli.verbose_ > 0) cout << "Create VBR from reordered CSR" << endl;
    vbmat.fill_from_CSR(cmat_A, partition, cli.block_size_);
    if (cli.verbose_ > 1) vbmat.print();


    intT A_rows = cmat_A.rows;
    intT A_cols = cmat_A.cols;
    intT B_cols = A_rows;
    intT B_rows = 100;
    intT C_cols = A_cols;
    intT C_rows = B_rows;

    //dense operand
    DataT* mat_B = new DataT[B_rows * B_cols];

    //TODO initialize B with random numbers, not 1s
    for (intT i = 0; i < B_rows * B_cols; i++)
    {
        mat_B[i] = 1.;
    }

    //******************************************
    //****VBR by dense MULTIPLICATION PHASE***
    //******************************************


    DataT_C* mat_C_VBR = new DataT_C[C_rows * C_cols]{ 0 }; //will store result of dense-dense multiplication

    //run the VBR-dense multiplications

    for (int i = -cli.warmup_; i < cli.exp_repetitions_; i++)
    {
        cublas_blockmat_multiply(vbmat, mat_B, B_cols, B_rows, mat_C_VBR, C_rows, dt, 8)
        //only saves non-warmup runs
        if (i >= 0) algo_times.push_back(dt);
    }

    //mean_time = mean(algo_times);
    //std_time = std_dev(algo_times);
    //algo_times.clear();

    //delete[] mat_A_gemm;



    //******************************************
    //****Sparse by dense MULTIPLICATION PHASE***
    //******************************************


    //TODO



    //******************************************
    //****VBR by dense MULTIPLICATION PHASE***
    //******************************************


    //TODO



}