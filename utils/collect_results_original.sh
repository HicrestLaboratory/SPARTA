#!/bin/bash

masks=( 1 16 64 256 )
taus=("01" "03" "06" "09")
centroids=(16 32 64 128)
col_reorder="0" #0: 1dim, 1: 2dim


SCRIPT_DIR=$(dirname "$0")
source "$SCRIPT_DIR/parse_args.sh" "$@"

mkdir -p "$output_dir"
output_file="${output_dir}/original_scramble${scramble}_bsize${block_size}.txt"
echo ${output_file}

> ${output_file}
echo "matrix_name rows cols nnz block_size scramble algo VBR_nzcount VBR_nzblocks_count VBR_average_height VBR_longest_row" >> ${output_file}


for matrix_folder in $(find "$matrix_dir" -mindepth 1 -maxdepth 1 -type d); do
    matrix_file=${matrix_folder}/$(basename "$matrix_folder").mtx

    echo " "
    echo "************************************************"
    matrix_name=$(basename "$matrix_file" .mtx)
    echo " Processing matrix: $matrix_name"

    #read matrix shape from mtx file
    header_line=$(grep -v '^%' "$matrix_file" | head -n 1)
    read -r rows cols nnz <<< "$header_line"

    matrix_result_dir="$result_dir/${matrix_name}"
    
    if [ ! -d "$matrix_result_dir" ]; then
        echo "Result folder not found: $matrix_result_dir"
        continue
    fi 
                
    if [ "${scramble}" == "1" ]; then
        echo "Scramble not implemented"
    continue #TODO use scramble function
    fi    
    

    params=$( echo "$matrix_name" "$rows" "$cols" "$nnz" "$block_size" "${scramble}" "original")
    club_outputs=$(./programs/general/Matrix_Analysis "$matrix_file" "$block_size")
    if [ $? -eq 0 ]; then
        # Only print and log the output if the command succeeds
        echo "$params $club_outputs" >> "${output_file}"
    else
        echo "Matrix_Analysis command failed for $matrix_file with block size $block_size"
    fi
echo "___________DONE"
done
