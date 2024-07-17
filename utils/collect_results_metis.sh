#!/bin/bash

metis_parts=( 1 2 4 8 16 32 64 128)
metis_objs=("edge-cut" "volume")

SCRIPT_DIR=$(dirname "$0")
source "$SCRIPT_DIR/parse_args.sh" "$@"

mkdir -p "$output_dir"
output_file=${output_dir}/metis_scramble${scramble}_bsize${block_size}.txt

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

    matrix_result_dir="$result_dir/${matrix_name}/GP"
    
    if [ ! -d "$matrix_result_dir" ]; then
        echo "Result folder not found: $matrix_result_dir"
        continue
    fi 


    for metis_part in "${metis_parts[@]}"; do

        matrix_result_dir="$result_dir/${matrix_name}/GP/${metis_part}"

        #check dir exists
        if [ ! -d "$matrix_result_dir" ]; then
        echo "Result folder not found for $matrix_name , $metis_parts parts"
        continue
        fi 
            
        for metis_obj in "${metis_objs[@]}"; do

            if [ "${scramble}" == "0" ]; then
                grouping_file="${matrix_result_dir}/${matrix_name}_metis_${metis_obj}_part${metis_part}.txt"
            else    
                grouping_file="${matrix_result_dir}/${matrix_name}-scramble_metis_${metis_obj}_part${metis_part}.txt"
            #check file exists
            fi

            if [ ! -f "$grouping_file" ]; then
            echo "Result file not found for $grouping_file"
            continue
            fi 

            params=$( echo "$matrix_name" "$rows" "$cols" "$nnz" "$block_size" "$scramble" "metis" "$metis_obj" "$metis_part")
            metis_outputs=$(./programs/general/Matrix_Analysis "$matrix_file" "$block_size" "$grouping_file" 1)
            if [ $? -eq 0 ]; then
                # Only print and log the output if the command succeeds
                echo "$params $metis_outputs" >> "${output_file}"
            else
                echo "Matrix_Analysis command failed for $matrix_file with block size $block_size"
            fi
        done
    done
echo "___________DONE"
done

