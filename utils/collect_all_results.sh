#!/bin/bash

# Default values for the variables
default_matrix_dir="../Clubs_tmp/matrices/toy"
default_result_dir="../../../Downloads/Club_result/Club_result/toy"
default_output_dir="results/results_2024"
default_block_size=64

# Function to display usage information
usage() {
    echo "Usage: $0 [-s script] [-m matrix_dir] [-r result_dir] [-o output_dir] [-b block_size]"
    echo "  -s script       Specify the script to run: clubs, metis, saad, all"
    echo "  -m matrix_dir   Specify the matrix directory (default: $default_matrix_dir)"
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
while getopts ":s:m:r:o:b:h" opt; do
    case "${opt}" in
        s)
            script_choice=${OPTARG}
            ;;
        m)
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

# Validate script_choice
if [[ "$script_choice" != "clubs" && "$script_choice" != "gp" && "$script_choice" != "all" && "$script_choice" != "saad" ]]; then
    echo "Invalid script choice: $script_choice. Allowed: clubs, metis, saad, all"
    usage
fi

mkdir -p "$output_dir"

# Array of scrambles
scrambles=(0 1)

echo "Matrix directory: $matrix_dir"
echo "Result directory: $result_dir"
echo "Output directory: $output_dir"

scrambles=(0 1)

# Run the scripts as executable commands
for scramble in "${scrambles[@]}"; do
    if [[ "$script_choice" == "all" || "$script_choice" == "clubs" ]]; then
        "./utils/collect_results_clubs.sh" -o ${output_dir} -b 64 -s "$scramble"
    fi
    
    if [[ "$script_choice" == "all" || "$script_choice" == "metis" ]]; then
        "./utils/collect_results_GP.sh" -o ${output_dir} -b 64 -s "$scramble"
    fi

done