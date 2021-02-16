#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <unistd.h> //getopt, optarg
#include <ctime>
#include <cmath>
#include <algorithm> //min_element, max_element




#include <vector>
#include <string>


#include "sparse_utilities.h"
#include "nvbs_utilities.h"
#include "comp_mats.h"

int main(int argc, char* argv[]) {

    int A_rows = 10;
    int A_cols = 20;
    int B_cols = 6;
    int block_size = 6;
    float mat_density = 0.5f;
    float row_density = 0.3f;
    float A_sparsity = 0.0f;
    intT block_cols = std::ceil((1.f*A_cols)/block_size);
    intT* col_part = new intT[block_cols + 1];
    partition(col_part, 0, A_cols, block_size);



    DataT* mat_A = new DataT[A_rows * A_cols]{ 0 };

    random_mat(mat_A, A_rows, A_cols, A_sparsity);
    
    std::cout << "\n A" << std::endl;
    matprint(mat_A, A_rows, A_cols, A_cols, 0);


    std::cout << "\n cmat A" << std::endl;
    CSR cmatA;
    convert_to_CSR(mat_A, A_rows, A_cols, 0, cmatA, 0);
    matprint(cmatA);

    std::cout << "\n cmat to vbmat" << std::endl;
    ncVBS vbmat_csr;
    convert_to_ncVBS(cmatA, vbmat_csr, block_cols, col_part);
    matprint(vbmat_csr);

    std::cout << "\n EQUALITY CHECK: vbmat and mat: " << equal(vbmat_csr, mat_A, A_rows, A_cols, A_cols, 0) << std::endl;

    /*
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
    */

    cleanVBS(vbmat_csr);
}