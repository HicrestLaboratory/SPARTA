#!/bin/bash

#************************************
#*************  INPUT  **************
#************************************************************************************************************************************************

# Default values
date_string=$(date +%d-%m-%Y)
experiment_folder="../../../Downloads/outputs_2024-08-01/SbatchMan/outputs/marzola"
collected_data_folder="results/results_$date_string/mult_data"
csv_folder="results/results_$date_string/mult_csv"

# Function to display usage
usage() {
    echo "Usage: $0 [-r experiment_dir] [-s storage_dir] [-c csv_dir]"
    echo "  -r  experiment_folder (default: $experiment_folder)"
    echo "  -o  collected_data_folder (default: $collected_data_folder)"
    echo "  -c  csv_folder (default: $csv_folder)"
    exit 1
}

# Parse options
while getopts ":r:s:c:o:" opt; do
    case ${opt} in
        r )
            experiment_folder=$OPTARG
            ;;
        o )
            collected_data_folder=$OPTARG
            ;;
        c )
            csv_folder=$OPTARG
            ;;
        \? )
            echo "Invalid option: -$OPTARG" 1>&2
            usage
            ;;
        : )
            echo "Invalid option: -$OPTARG requires an argument" 1>&2
            usage
            ;;
    esac
done
shift $((OPTIND -1))

# Print the values to check
echo "Experiment Folder: $experiment_folder"
echo "Data collection folder: $collected_data_folder"
echo "CSV Folder: $csv_folder"

routines=("spmmcsr" "spmvcsr" "spmmbsr" "spmvbsr")
methods=("clubs" "metis-edge-cut" "metis-volume" "saad" "original" "denseAMP")
mkdir -p "$collected_data_folder"
mkdir -p "$csv_folder"



for routine in "${routines[@]}";do
    echo "**********   COLLECTING MULTIPLICATION DATA FOR $routine"
    ./utils/rearrange_routine_mult_results.sh -x -r "$experiment_folder" -o "$collected_data_folder" -u "$routine"
    echo "**********   CREATING MULTIPLICATION CSV FOR $routine $method"
    ./utils/collect_routine_mult_results.sh -r "$collected_data_folder" -c "$csv_folder" -u "$routine" -m "ALL"
done


