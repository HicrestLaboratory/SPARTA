#pragma once
#include<vector>

typedef long int intT;
typedef float DataT; //precision for input matrices entries
typedef float DataT_C; //precision for result matrix entries
//types for the two multiplied matrices and the result matrix. available choices are: 
//          INT_8, INT_32
//          FLOAT_32, FLOAT_32

typedef float (*distFunc)(intT*,intT,intT*,intT);
typedef float (*distFuncQuot)(intT*,intT,intT*,intT, intT );
typedef float (*distFuncGroup)(std::vector<intT>,intT,intT*,intT,intT,intT);

enum BlockingType {iterative, iterative_structured, fixed_size, iterative_clocked, iterative_queue, iterative_max_size};

enum MultiplicationAlgo {NO_MULT, cublas_gemm, cusparse_spmm, cusparse_bellpack, cublas_vbr, cublas_vbr_fixed};