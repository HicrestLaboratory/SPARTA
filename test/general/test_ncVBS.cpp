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

    ncVBS vbmat;
    int A_rows = 4;
    int A_cols = 3;
    int B_cols = 5;
    int block_size = 3;
    float mat_density = 0.5f;
    float row_density = 0.3f;
    float A_sparsity = 0.5f;

    DataT* mat_A = new DataT[A_rows * A_cols]{ 0 };
    
    random_mat(mat_A, A_rows, A_cols, A_sparsity);
    
    std::cout << "A" << std::endl;
    matprint(mat_A, A_rows, A_cols, A_cols, 0);


    std::cout << "VBMAT_A" << std::endl;
    block_cols = A_cols / block_size;
    intT* col_part = new intT[block_cols + 1];
    partition(col_part, 0, A_cols, block_size);
    convert_to_ncVBS(mat_A, A_rows, A_cols, 0, A_cols, vbmat, block_cols, col_part);
    matprint(vbmat);


    int B_rows = A_cols;
    DataT* mat_B = new DataT[B_rows * B_cols]{ 0 };

    random_mat(mat_B, B_rows, B_cols, A_sparsity);

    std::cout << "B" << std::endl;
    matprint(mat_B, B_rows, B_cols, B_cols, 0);

    int C_rows = A_rows;
    int C_cols = B_cols;
    DataT* mat_C = new DataT[C_rows * C_cols]{ 0 };
    multiply(mat_A, A_rows, A_cols, 0, A_cols, mat_B, B_cols, 0, B_cols, mat_C, C_cols, 0);

    std::cout << "C" << std::endl;
    matprint(mat_C, C_rows, C_cols, C_cols, 0);


}