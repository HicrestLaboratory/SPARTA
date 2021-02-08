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
    int mat_rows = 10;
    int mat_cols = 18;
    int block_size = 3;
    float mat_density = 0.5f;
    float row_density = 0.3f;

    random_ncVBS(vbmat, mat_rows, mat_cols, block_size, mat_density, row_density);
    matprint(vbmat);

    DataT* mat = new DataT[mat_rows * mat_cols]{ 0 };
    convert_to_mat(vbmat, mat, 0);
    matprint(mat, mat_rows, mat_cols, mat_cols, 0);

    std::cout << "MAT2" << std::endl;
    DataT* mat2 = new DataT[mat_rows * mat_cols]{ 0 };
    mat2[1] = 1;
    mat2[3] = 1;
    mat2[10] = 1;
    mat2[20] = 1;
    mat2[25] = 1;
    mat2[26] = 1;
    mat2[31] = 1;
    mat2[100] = 1;
    mat2[123] = 1;


    ncVBS vbmat2;
    matprint(mat2, mat_rows, mat_cols, mat_cols, 0);

    int block_cols = mat_cols / block_size;
    int* col_part = new int[block_cols];
    partition(col_part, 0, mat_cols, block_size);

    
    std::cout << "VBMAT2" << std::endl;

    convert_to_ncVBS(mat2, mat_rows, mat_cols, 0, mat_cols, vbmat2, block_cols, col_part);
    matprint(vbmat2);

    DataT* mat3 = new DataT[mat_rows * mat_cols]{ 0 };
    convert_to_mat(vbmat2, mat3, 0);


    std::cout << "MAT3" << std::endl;

    matprint(mat3, mat_rows, mat_cols, mat_cols, 0);
    
    std::cout << "EQUAL?" << equal(mat_rows, mat_cols, mat2, mat_cols, 0, mat3, mat_cols, 0, 0.1f) << std::endl;;

}