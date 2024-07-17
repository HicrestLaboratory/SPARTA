#!/bin/bash

# Function to display help message
usage() {
    echo "Usage: $0 [-f matrix_dir] [-r result_dir] [-o output_dir] [-b block_size] [-x metis_obj] [-p metis_part] [-t tau] [-m mask] [-c centroids] [-s scramble]"
    exit 1
}

# Default values for the variables
matrix_dir="../Clubs_tmp/matrices/medium"
result_dir="../../../Downloads/Club_result/Club_result/medium"
output_dir="results/results_2024/medium"
block_size=64
metis_obj="edge-cut"
metis_part=128
tau="01"
mask="64"
centroids="64"
scramble="0"

# Parse command-line arguments
while getopts ":f:r:b:o:x:p:t:m:c:s:h" opt; do
    case "${opt}" in
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
        x)
            metis_obj=${OPTARG}
            ;;
        p)
            metis_part=${OPTARG}
            ;;
        t)
            tau=${OPTARG}
            ;;
        m)
            mask=${OPTARG}
            ;;
        c)
            centroids=${OPTARG}
            ;;
        s)
            scramble=${OPTARG}
            ;;
        h)
            usage
            ;;
        \?)
            echo "Invalid option: -${OPTARG}" >&2
            usage
            ;;
        :)
            echo "Option -${OPTARG} requires an argument." >&2
            usage
            ;;
    esac
done
shift $((OPTIND -1))

# Export variables to be used by other scripts
export matrix_dir
export result_dir
export output_dir
export block_size
export metis_obj
export metis_part
export tau
export centroids
export mask
export scramble