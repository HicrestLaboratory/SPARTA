#!/bin/bash

# Default values for the variables
default_matrix_dir="../Clubs_tmp/matrices/toy"
default_result_dir="../../../Downloads/Club_result/Club_result/toy"
default_output_dir="results/results_2024/"
default_block_size=64
valid_choices=( "all" "clubs" "metis" "saad" "original" "patoh" )

# Function to display usage information
usage() {
    echo "Usage: $0 [-s script] [-m matrix_dir] [-r result_dir] [-o output_dir] [-b block_size]"
    echo "  -s script       Specify the script to run: ${valid_choices[*]}"
    echo "  -f matrix_dir   Specify the matrix directory (default: $default_matrix_dir)"
    echo "  -r result_dir   Specify the result directory (default: $default_result_dir)"
    echo "  -o output_dir   Specify the output directory (default: $default_output_dir)"
    echo "  -b block_size   Specify the block size (default: $block_size)"
    exit 1
}

# Initialize variables with default values
script_choice="all"
matrix_dir="$default_matrix_dir"
result_dir="$default_result_dir"
output_dir="$default_output_dir"
block_size="$default_block_size"

# Parse command-line options
while getopts ":s:f:r:o:b:h" opt; do
    case "${opt}" in
        s)
            script_choice=${OPTARG}
            ;;
        f)
            matrix_dir=${OPTARG}
            ;;
        r)
            result_dir=${OPTARG}
            ;;
        o)
            output_dir=${OPTARG}
            ;;
        b)
            block_size=${OPTARG}
            ;;
        h)
            usage
            ;;
        *)
            usage
            ;;
    esac
done

# Check if script_choice is valid
if [[ ! " ${valid_choices[@]} " =~ " ${script_choice} " ]]; then
    echo "Invalid script choice: '$script_choice'. Allowed choices are: ${valid_choices[*]}"
    usage
fi

mkdir -p "$output_dir"

# Array of scrambles

echo "Matrix directory: $matrix_dir"
echo "Result directory: $result_dir"
echo "Output directory: $output_dir"

args=(-f "$matrix_dir" -r "$result_dir" -o "$output_dir" -b "$block_size")

# Run the scripts as executable commands
scramble="0"


if [[ "$script_choice" == "all" || "$script_choice" == "original" ]]; then

    echo "                        ***************************************"
    echo "************************************ ORIGINAL *****************"
    echo "                        ***************************************"

    "./utils/collect_results_original.sh" "${args[@]}" -s "$scramble"
fi


if [[ "$script_choice" == "all" || "$script_choice" == "clubs" ]]; then
    echo "                        ***************************************"
    echo "*************************************** CLUBS *****************"
    echo "                        ***************************************"
    "./utils/collect_results_clubs.sh" "${args[@]}" -s "$scramble"

fi


if [[ "$script_choice" == "all" || "$script_choice" == "metis" ]]; then
    echo "                        ***************************************"
    echo "*************************************** METIS *****************"
    echo "                        ***************************************"
    "./utils/collect_results_metis.sh" "${args[@]}" -s "$scramble"
fi


if [[ "$script_choice" == "all" || "$script_choice" == "denseamp" ]]; then
    echo "                        ***************************************"
    echo "*************************************** DENSEAMP *****************"
    echo "                        ***************************************"
    "./utils/collect_results_denseamp.sh" "${args[@]}" -s "$scramble"
fi


if [[ "$script_choice" == "all" || "$script_choice" == "patoh" ]]; then
    echo "                        ***************************************"
    echo "*************************************** PATOH *****************"
    echo "                        ***************************************"
    "./utils/collect_results_patoh.sh" "${args[@]}" -s "$scramble"
fi


if [[ "$script_choice" == "all" || "$script_choice" == "saad" ]]; then
    echo "                        ***************************************"
    echo "*************************************** SAAD *****************"
    echo "                        ***************************************"
    "./utils/collect_results_saad.sh" "${args[@]}" -s "$scramble"
fi