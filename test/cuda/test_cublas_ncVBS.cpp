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
#include "cuda_utilities_ncVBS.h"

using namespace std;


int main(int argc, char* argv[]) {

	ncVBS vbmat;
	int A_rows = 100;
	int A_cols = 20;
	int B_cols = 20;
	int block_size = 5;
//	float mat_density = 0.1f;
//	float row_density = 0.1f;
	float A_sparsity = 0.1f;
	float B_sparsity = 0.9f;
	int streams = 32;


    opterr = 0;
    char c;
    while ((c = getopt(argc, argv, "a:b:i:q:e:f:m:n:p:r:k:s:v:w:S:")) != -1)
        switch (c)
        {
		case 'm':
			A_rows = stoi(optarg);
			break;
		case 'n':
			B_cols = stoi(optarg);
			break;
		case 'k':
			A_cols = stoi(optarg);
			break;
		case 'p':
			block_size = stoi(optarg);
			break;
		case 's':
			streams = stoi(optarg);
        case '?':
            fprintf(stderr, "Option -%c does not exists, or requires an argument.\n", optopt);
            return 1;
        default:
            abort();
        }


	int C_rows = A_rows;
	int C_cols = B_cols;

	//CREATING A RANDOM MATRIX (A)
	DataT* mat_A = new DataT[A_rows * A_cols]{ 0 };
	random_mat(mat_A, A_rows, A_cols, A_sparsity);
	std::cout << "\n A" << std::endl;
	//matprint(mat_A, A_rows, A_cols, A_cols, 0);

	//CONVERTING TO VBMAT
	std::cout << "\n VBMAT_A" << std::endl;
	intT block_cols = A_cols / block_size;
	intT* col_part = new intT[block_cols + 1];
	partition(col_part, 0, A_cols, block_size);
	convert_to_ncVBS(mat_A, A_rows, A_cols, 0, A_cols, vbmat, block_cols, col_part);
	//matprint(vbmat);

	//CREATING A RANDOM MATRIX (B)
	int B_rows = A_cols;
	DataT* mat_B = new DataT[B_rows * B_cols]{ 0 };
	random_mat(mat_B, B_rows, B_cols, B_sparsity);
	std::cout << "\n B created" << std::endl;
	//matprint(mat_B, B_rows, B_cols, B_cols, 0);


	bool do_serial = false;

	DataT* mat_C = new DataT[C_rows * C_cols]{ 0 };
	if (do_serial)
	{
		//MULTIPLYING SERIAL
		multiply(mat_A, A_rows, A_cols, 0, A_cols, mat_B, B_cols, 0, B_cols, mat_C, C_cols, 0);
		std::cout << "\n C" << std::endl;
		//matprint(mat_C, C_rows, C_cols, C_cols, 0);


		//MULTIPLYING WITH VBMAT
		std::cout << "\n multiplying with vbmat" << std::endl;
		DataT* mat_C2 = new DataT[C_rows * C_cols]{ 0 };
		multiply(vbmat, mat_B, B_cols, 0, B_cols, mat_C2, C_cols, 0);
		//matprint(mat_C2, C_rows, C_cols, C_cols, 0);
		std::cout << "EQUALITY CHECK: MULTIPLICATION: " << equal(C_rows, C_cols, mat_C, C_cols, 0, mat_C2, C_cols, 0, 0.00001f) << std::endl;
	}

	//MULTIPLYING WITH VBMAT CUBLAS
	std::cout << "\n multiplying with vbmat CUBLAS" << std::endl;
	DataT* mat_C_vbs_cublas = new DataT[C_rows * C_cols]{ 0 };
	float* dt = new float[6];
	cublas_ncVBS_multiply(vbmat, mat_B, B_cols, B_cols, mat_C_vbs_cublas, C_cols, dt, streams, streams);
	//matprint(mat_C3, C_rows, C_cols, C_cols, 0);

	if (do_serial) std::cout << "EQUALITY CHECK: MULTIPLICATION: " << equal(C_rows, C_cols, mat_C, C_cols, 0, mat_C_vbs_cublas, C_cols, 0, 0.00001f) << std::endl;
	std::cout << "TIME MEASUREMENTS:" << std::endl;
	for (int i = 0; i < 6; i++) std::cout << dt[i] << " ";
	std::cout << std::endl;

	//MULTIPLYING WITH CUSPARSE

	//TODO do not start with array but create directly the CSR?

	CSR cmat_A;
	convert_to_CSR(mat_A, A_rows, A_cols, 0, cmat_A, 0);
	//prepare the cusparse CSR format
	int nnz = count_nnz(cmat_A);
	int* csrRowPtr = new int[cmat_A.rows + 1];
	int* csrColInd = new int[nnz];
	float* csrVal = new float[nnz];
	prepare_cusparse_CSR(cmat_A, csrRowPtr, csrColInd, csrVal);

	float cusparse_dt;
	DataT* mat_C_cusparse = new DataT[C_rows * C_cols]{ 0 };
	std::cout << "\n multiplying with CUSPARSE" << std::endl;
	cusparse_gemm_custom(cmat_A.rows, cmat_A.cols, nnz, csrRowPtr, csrColInd, csrVal, mat_B, B_cols, B_rows, mat_C_cusparse, C_rows, 1.0f, 0.0f, cusparse_dt);
	if (do_serial) std::cout << "EQUALITY CHECK: MULTIPLICATION: " << equal(C_rows, C_cols, mat_C, C_cols, 0, mat_C_cusparse, C_cols, 0, 0.00001f) << std::endl;
	std::cout << "TIME MEASUREMENT:" << cusparse_dt << std::endl;

	cleanVBS(vbmat);
	cleanCSR(cmat_A);
}
 
 