#!/bin/bash

# PaToH partition sizes
parts=( 1 2 4 8 16 32 64 128 )

SCRIPT_DIR=$(dirname "$0")
source "$SCRIPT_DIR/parse_args.sh" "$@"

mkdir -p "$output_dir"
output_file="${output_dir}/patoh_scramble${scramble}_bsize${block_size}.txt"

> "$output_file"
echo "matrix_name rows cols nnz block_size scramble algo patoh_part VBR_nzcount VBR_nzblocks_count VBR_average_height VBR_longest_row" >> "$output_file"

for matrix_folder in $(find "$matrix_dir" -mindepth 1 -maxdepth 1 -type d); do
    matrix_file="${matrix_folder}/$(basename "$matrix_folder").mtx"

    echo "************************************************"
    matrix_name=$(basename "$matrix_file" .mtx)
    echo " Processing matrix: $matrix_name"
    echo "************************************************"

    # Read matrix shape from mtx file
    header_line=$(grep -v '^%' "$matrix_file" | head -n 1)
    read -r rows cols nnz <<< "$header_line"

    matrix_result_dir="$result_dir/${matrix_name}/patoh"
    
    if [ ! -d "$matrix_result_dir" ]; then
        echo "Result folder not found: $matrix_result_dir"
        continue
    fi 

    for patoh_part in "${parts[@]}"; do

        # Correct directory for each partition
        matrix_part_dir="$matrix_result_dir/${patoh_part}"

        # Check if directory exists
        if [ ! -d "$matrix_part_dir" ]; then
            echo "Result folder not found for $matrix_name, ${patoh_part} parts"
            continue
        fi 

        grouping_file="${matrix_part_dir}/${matrix_name}_cutpart_k${patoh_part}_default_s1_partvec.txt"

        # Check if file exists
        if [ ! -f "$grouping_file" ]; then
            echo "Result file not found: $grouping_file"
            continue
        fi 

        # Correctly use $patoh_part in params
        params=$(echo "$matrix_name" "$rows" "$cols" "$nnz" "$block_size" "$scramble" "patoh" "$patoh_part")
        patoh_outputs=$(./programs/general/Matrix_Analysis "$matrix_file" "$block_size" "$grouping_file" 1)
        
        if [ $? -eq 0 ]; then
            # Only print and log the output if the command succeeds
            echo "$params $patoh_outputs" >> "$output_file"
        else
            echo "Matrix_Analysis command failed for $matrix_file with block size $block_size"
        fi

    done
done

echo "___________DONE"
