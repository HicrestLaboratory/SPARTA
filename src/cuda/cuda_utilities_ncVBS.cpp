//CUDA Utilities and system includes
#include <assert.h>
#include <helper_string.h>  // helper for shared functions common to CUDA Samples

// CUDA runtime
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>

// CUDA and CUBLAS functions
#include <helper_functions.h>
#include <helper_cuda.h>

//others
#include<stdio.h>
#include <iostream>
#include <typeinfo>

//includes from SPARTA
#include "cuda_utilities.h"
#include "cuda_utilities_ncVBS.h"
#include "comp_mats.h"
#include "sparse_utilities.h"

/*ONLY WORKS WITH 
    DataT = float, double, int;
    intT = int
    unsigned_intT = unsigned int;
*/



void cublas_ncVBS_multiply(const ncVBS& vbmatA, DataT* B, int B_cols, int B_lead_dim, DataT* C, int C_lead_dim, float& dt)
{
    //NOT IMPLEMENTED YET
}