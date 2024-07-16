#!/bin/bash

# Paths (you can adjust these as needed)
matrix_dir="../Clubs_tmp/matrices/toy"
result_dir="../../../Downloads/Club_result/Club_result/toy"
block_size=100
metis_obj="edge-cut"
metis_part=128

# Loop through each matrix file in the matrix directory
for matrix_file in $(find "$matrix_dir" -type f -name "*.mtx"); do
    echo "************************************************"
    echo " Processing file: $matrix_file"
    matrix_name=$(basename "$matrix_file" .mtx)
    
    matrix_result_dir="$result_dir/${matrix_name}/GP/${metis_part}"
    echo "Searching results in ${matrix_result_dir}"
    if [ ! -d "$matrix_result_dir" ]; then
	echo "Result folder not found for $matrix_name"
	continue
    fi 
    result_file="${matrix_result_dir}/${matrix_name}_metis_${metis_obj}_part${metis_part}.txt"
    ./programs/general/Matrix_Analysis "$matrix_file" "$block_size" "$result_file" 1
done
