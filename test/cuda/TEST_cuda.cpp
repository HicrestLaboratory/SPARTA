#include "cuda_utilities.h"
#include "matrices.h"
#include "blocking.h"
#include "utilies.h"

//keep track of time
float dt;
vec_d algo_times;
float mean_time;
float std_time;

//import CSR
CLineReader cli(argc, argv);
if (cli.verbose_ > 0) cli.print();
CSR cmat_A(cli); //sparse operand
BlockingEngine bEngine(cli);

//dense operand
DataT* mat_B = new DataT[B_rows * B_cols];

//TODO initialize B with random numbers, not 1s
for (intT i = 0; i < B_rows * B_cols; i++)
{
    mat_B[i] = 1.;
}

//******************************************
//****Dense by dense MULTIPLICATION PHASE***
//******************************************


//TODO convert csr cmat_A to dense mat_A_gemm
DataT* mat_A_gemm = new DataT[A_rows * A_cols]{ 0 };

DataT_C* mat_Cgemm = new DataT_C[C_rows * C_cols]{ 0 }; //will store result of dense-dense multiplication

//run the dense-dense multiplications TODO

//for (int i = -cli.warmup_; i < cli.exp_repetitions_; i++)
//{
//    cublas_gemm_custom(mat_A_gemm, A_rows, A_cols, A_rows, mat_B, B_cols, B_rows, mat_Cgemm, C_rows, 1, 0, dt);
//    //only saves non-warmup runs
//    if (i >= 0) algo_times.push_back(dt);
//}

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