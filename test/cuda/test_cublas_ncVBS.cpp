#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <string>
#include <unistd.h>
#include <math.h>
#include <typeinfo>
#include <iterator>
#include <algorithm>

// Utilities and system include
#include <assert.h>
#include <helper_string.h>  // helper for shared functions common to CUDA Samples

// CUDA runtime
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>


// CUDA and CUBLAS functions
#include <helper_functions.h>
#include <helper_cuda.h>

#include "cuda_utilities.h"
#include "sparse_utilities.h"
#include "nvbs_utilities.h"
#include "comp_mats.h"



int main(int argc, char* argv[]) {

	ncVBS vbmat;
	int A_rows = 4;
	int A_cols = 9;
	int B_cols = 5;
	int block_size = 3;
	float mat_density = 0.5f;
	float row_density = 0.3f;
	float A_sparsity = 0.5f;

	DataT* mat_A = new DataT[A_rows * A_cols]{ 0 };

	random_mat(mat_A, A_rows, A_cols, A_sparsity);

	std::cout << "\n A" << std::endl;
	matprint(mat_A, A_rows, A_cols, A_cols, 0);


	std::cout << "\n VBMAT_A" << std::endl;
	intT block_cols = A_cols / block_size;
	intT* col_part = new intT[block_cols + 1];
	partition(col_part, 0, A_cols, block_size);
	convert_to_ncVBS(mat_A, A_rows, A_cols, 0, A_cols, vbmat, block_cols, col_part);
	matprint(vbmat);
	std::cout << "\n EQUALITY CHECK: vbmat and mat: " << equal(vbmat, mat_A, A_rows, A_cols, A_cols, 0) << std::endl;


	int B_rows = A_cols;
	DataT* mat_B = new DataT[B_rows * B_cols]{ 0 };

	random_mat(mat_B, B_rows, B_cols, A_sparsity);

	std::cout << "\n B" << std::endl;
	matprint(mat_B, B_rows, B_cols, B_cols, 0);

	int C_rows = A_rows;
	int C_cols = B_cols;
	DataT* mat_C = new DataT[C_rows * C_cols]{ 0 };
	multiply(mat_A, A_rows, A_cols, 0, A_cols, mat_B, B_cols, 0, B_cols, mat_C, C_cols, 0);

	std::cout << "\n C" << std::endl;
	matprint(mat_C, C_rows, C_cols, C_cols, 0);

	std::cout << "\n C multiplied with vbmat" << std::endl;
	DataT* mat_C2 = new DataT[C_rows * C_cols]{ 0 };
	multiply(vbmat, mat_B, B_cols, 0, B_cols, mat_C2, C_cols, 0);
	matprint(mat_C2, C_rows, C_cols, C_cols, 0);

	std::cout << "\n EQUALITY CHECK: MULTIPLICATION: " << equal(C_rows, C_cols, mat_C, C_cols, 0, mat_C2, C_cols, 0, 0.00001f) << std::endl;


	DataT* mat_C3 = new DataT[C_rows * C_cols]{ 0 };
	float dt = 0;
	cublas_ncVBS_multiply(vbmat, mat_B, B_cols, B_cols, mat_C3, C_cols, dt);



}
 
 