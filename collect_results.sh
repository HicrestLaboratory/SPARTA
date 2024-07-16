#!/bin/bash

# Paths (you can adjust these as needed)
matrix_dir="../Clubs/matrices"
result_dir="../../../Downloads/Club_result/Club_result"


# Loop through each matrix file in the matrix directory
for matrix_file in "$matrix_dir"/*.mtx; do
    matrix_name=$(basename "$matrix_file" .mtx)
    
    matrix_result_dir="$result_dir/$matrix_name"
    if [! -d "$matrix_result_dir" ]; then
	echo "Result folder not found for $matrix_name"
	continue
        # Loop through each technique directory within the result directory
    technique_folder_dir=""
        # Loop through each file in the technique directory
        for result_file in "$technique_dir"/*; do
            # Run the program with the given parameters
            ./programs/general/Matrix_Analysis "$matrix_file" 105 "$result_file" 1
        done
    done
done

