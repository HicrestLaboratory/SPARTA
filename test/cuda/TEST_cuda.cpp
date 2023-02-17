#include "cuda_utilities.h"
#include "matrices.h"
#include "blocking.h"
#include "utilities.h"
#include <iostream>
#include <random>
#include <algorithm> //fill, equal


using namespace std;



int main(int argc, char* argv[])
{
    //keep track of time
    float dt;
    vector<double> algo_times;
    float mean_time;
    float std_time;

    bool check_results = true;


    //import CSR
    CLineReader cli(argc, argv);
    if (cli.verbose_ > 0) cli.print();
    CSR cmat_A(cli); //sparse operand
    BlockingEngine bEngine(cli);
    //bEngine.blocking_algo = fixed_size; //set blocking algo to fixed_size
    //bEngine.row_block_size = bEngine.col_block_size;


    VBR vbmat; 
    if (cli.verbose_ > 0) cout << "Blocking the matrix" << endl;
    bEngine.GetGrouping(cmat_A);

    if (cli.verbose_ > 0) cout << "Creating VBR" << endl;
    vbmat.fill_from_CSR_inplace(cmat_A, bEngine.grouping_result, cli.col_block_size_);
    if (cli.verbose_ > 1) vbmat.print();

    intT A_rows = cmat_A.rows;
    intT A_cols = cmat_A.cols;
    intT B_rows = A_cols;
    intT B_cols = 128;
    intT C_cols = A_cols;
    intT C_rows = B_rows;

    //dense operand
    DataT* mat_B = new DataT[B_rows * B_cols];

    random_device rd;

    mt19937 e2(rd());
    uniform_real_distribution<> dist(0, 1);

    for (int n = 0; n < B_rows*B_cols; n++) 
    {
        mat_B[n] = 1.;//dist(e2);
    }
    if (cli.verbose_ > 1) 
    {
        cout << "Matrix B: " << B_rows << " X " << B_cols << endl;
        print_mat(mat_B, B_rows, B_cols, B_cols);
    }


    //******************************************
    //*****SERIAL CSR MULTIPLICATION PHASE******
    //******************************************

    DataT* mat_C_serial = new DataT[C_rows*C_cols];

    if (cli.verbose_ > 0) cout << "CSR Serial Multiplication" << endl;
    cmat_A.multiply(mat_B, B_cols, mat_C_serial);
    if (cli.verbose_ > 1) 
    {
        cout << "Serial CSR: multiplication result" << endl;
        print_mat(mat_C_serial, C_rows, C_cols, C_cols);
    }


    //******************************************
    //****VBR by dense MULTIPLICATION PHASE***
    //******************************************

    DataT_C* mat_C_VBR = new DataT_C[C_rows * C_cols]{ 0 }; //will store result of dense-dense multiplication
    algo_times.clear();
    //run the VBR-dense multiplications
    for (int i = -cli.warmup_; i < cli.exp_repetitions_; i++)
    {
        fill(mat_C_VBR, mat_C_VBR + C_cols*C_rows, 0);
        //cublas_blockmat_multiply(vbmat, mat_B, B_cols, B_rows, mat_C_VBR, C_rows, dt, 8);
        //only saves non-warmup runs
        if (i >= 0) algo_times.push_back(dt);
        bool equality_check = equal(mat_C_VBR, mat_C_VBR + C_cols*C_rows, mat_C_serial);
        cout << "CORRECTNESS CHECK: " << equality_check << endl;
    }

    cout << "TIME: " << avg(algo_times) << endl;

    //mean_time = mean(algo_times);
    //std_time = std_dev(algo_times);
    //algo_times.clear();

    delete[] mat_C_VBR;



    //******************************************
    //****Sparse by dense MULTIPLICATION PHASE***
    //******************************************


    //TODO



    //******************************************
    //****VBR by dense MULTIPLICATION PHASE***
    //******************************************


    //TODO



}