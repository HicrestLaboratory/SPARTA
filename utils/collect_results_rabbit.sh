#!/bin/bash
SCRIPT_DIR=$(dirname "$0")
source "$SCRIPT_DIR/parse_args.sh" "$@"

mkdir -p "$output_dir"
output_file=${output_dir}/rabbit_bsize${block_size}.txt

> "$output_file"
echo "matrix_name rows cols nnz block_size scramble algo metis_obj metis_part VBR_nzcount VBR_nzblocks_count VBR_average_height VBR_longest_row" >> ${output_file}


for matrix_folder in $(find "$matrix_dir" -mindepth 1 -maxdepth 1 -type d); do
    matrix_file=${matrix_folder}/$(basename "$matrix_folder").mtx

    echo "************************************************"
    matrix_name=$(basename "$matrix_file" .mtx)
    echo " Processing matrix: $matrix_name"
    echo "************************************************"

    #read matrix shape from mtx file
    header_line=$(grep -v '^%' "$matrix_file" | head -n 1)
    read -r rows cols nnz <<< "$header_line"

    matrix_result_dir="$result_dir/${matrix_name}/rabbit"
    
    if [ ! -d "$matrix_result_dir" ]; then
        echo "Result folder not found: $matrix_result_dir"
        continue
    fi 
        
    grouping_file="${matrix_result_dir}/${matrix_name}.g"
    
    if [ ! -f "$grouping_file" ]; then
        echo "Result file not found for $grouping_file"
        continue
    fi 

    params=$( echo "$matrix_name" "$rows" "$cols" "$nnz" "$block_size" "$scramble" "rabbit")
    metis_outputs=$(./programs/general/Matrix_Analysis "$matrix_file" "$block_size" "$grouping_file" 1)
    if [ $? -eq 0 ]; then
        # Only print and log the output if the command succeeds
        echo "$params $metis_outputs" >> "${output_file}"
    else
        echo "Matrix_Analysis command failed for $matrix_file with block size $block_size"
    fi
echo "___________DONE"
done

