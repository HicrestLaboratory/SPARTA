#!/bin/bash

taus=("01" "03" "06" "09")
col_reorder="0" #0: 1dim, 1: 2dim

echo "ARGS: " "$@"
SCRIPT_DIR=$(dirname "$0")
source "$SCRIPT_DIR/parse_args.sh" "$@"

mkdir -p "$output_dir"
output_file="${output_dir}/saad_scramble${scramble}_bsize${block_size}.txt"
echo ${output_file}

> ${output_file}
echo "matrix_name rows cols nnz block_size scramble algo tau VBR_nzcount VBR_nzblocks_count VBR_average_height VBR_longest_row" >> ${output_file}


for matrix_folder in $(find "$matrix_dir" -mindepth 1 -maxdepth 1 -type d); do
    matrix_file=${matrix_folder}/$(basename "$matrix_folder").mtx

    echo "************************************************"
    matrix_name=$(basename "$matrix_file" .mtx)
    echo " Processing matrix: $matrix_name"
    echo "************************************************"

    #read matrix shape from mtx file
    header_line=$(grep -v '^%' "$matrix_file" | head -n 1)
    read -r rows cols nnz <<< "$header_line"

    matrix_result_dir="$result_dir/${matrix_name}"
    
    if [ ! -d "$matrix_result_dir" ]; then
        echo "Result folder not found: $matrix_result_dir"
        continue
    fi 

    for tau in "${taus[@]}"; do
     
        if [ "${scramble}" == "1" ]; then
        grouping_file="${matrix_result_dir}/${matrix_name}-scramble_tau${tau}.g"
        else    
        grouping_file="${matrix_result_dir}/${matrix_name}_tau${tau}.g"
        #check file exists
        fi
        if [ ! -f "$grouping_file" ]; then
        echo "Result file not found for $grouping_file"
        continue
        fi 

        params=$( echo "$matrix_name" "$rows" "$cols" "$nnz" "$block_size" "${scramble}" "saad" "$tau")
        saad_outputs=$(./programs/general/Matrix_Analysis "$matrix_file" "$block_size" "$grouping_file" "$col_reorder")
        if [ $? -eq 0 ]; then
            # Only print and log the output if the command succeeds
            echo "$params $saad_outputs" >> "${output_file}"
        else
            echo "Matrix_Analysis command failed for $matrix_file with block size $block_size"
        fi
    done
echo "___________DONE"
done