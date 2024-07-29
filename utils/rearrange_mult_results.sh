#!/bin/bash

#************************************
#*************  INPUT  **************
#************************************************************************************************************************************************

# Default values
root_dir="../../../Downloads/outputs_2024-07-26/SbatchMan/outputs/marzola"
routine="spmmcsr"
output_dir="results/results_2024/mult_data"

# Function to display usage
usage() {
    echo "Usage: $0 [-r root_dir] [-u routine] [-m method]"
    echo "  -r  Root directory (default: $root_dir)"
    echo "  -u  Routine (default: $routine)"
    exit 1
}

# Parse options
while getopts ":r:u:" opt; do
    case ${opt} in
        r )
            root_dir=$OPTARG
            ;;
        u )
            routine=$OPTARG
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

# Define directories based on the parsed or default values
routine_dir="$root_dir/$routine/finished"
routine_gp_dir="$root_dir/${routine}gp/finished"

# Print the values to check
echo "Root directory: $root_dir"
echo "Routine: $routine"
echo "Routine directory: $routine_dir"
echo "Routine GP directory: $routine_gp_dir"

if [ ! -d "$routine_dir" ]; then
  echo "ERROR: COULD NOT FIND RESULT DIRECTORY $routine_dir"
  exit 1
fi

# Output CSV file
out_routine_dir=$output_dir/$routine
echo $out_routine_dir
mkdir -p "$out_routine_dir"

algos=( "metis-edge-cut" "metis-volume" "clubs" "denseAMP" "saad" "original")
for algo in "${algos[@]}"; do
  mkdir -p "$out_routine_dir/$algo"
done
#************************************
#************************************
#************************************************************************************************************************************************

#************************************************************************************************************************************************
#************* UTILS    *************
#************************************

# Function to display progress bar
show_progress() {
    local current=$1
    local total=$2
    local progress=$((current * 100 / total))
    local done=$((progress * 4 / 10))
    local left=$((40 - done))
    local fill=$(printf "%${done}s")
    local empty=$(printf "%${left}s")
    printf "\rProgress : [${fill// /#}${empty// /-}] $progress%%"
}

find_file_type() {
  local file_name=$1

  if [[ $file_name == *"msk"* && $file_name == *"cent"* && $file_name == *"tau"* ]]; then
    file_type="clubs"
  elif [[ $file_name == *"metis"* ]]; then 
    if [[ $file_name == *"edge"* ]];then
      file_type="metis-edge-cut"
    else
      file_type="metis-volume"
    fi
  elif [[ $file_name == *"Saad"* ]]; then
    file_type="saad"
  elif [[ $file_name == *"DenseAMP"* ]]; then
    file_type="denseAMP"
  else
    file_type="original"
  fi

  echo $file_type
}

process_file()
{
  local file_path=$1

  local file_name=$(basename "$file_path")

  local file_type=$(find_file_type "$file_name")

  local file_dir=$out_routine_dir/$file_type

  cp $file_path $file_dir

}


#************************************
#************************************
#************************************************************************************************************************************************

#************************************************************************************************************************************************

#-----------------------------------------



#find files
total_files=$(find "$routine_dir" -name '*.out' | wc -l)
gp_files=$(find "$routine_gp_dir" -name '*.out' | wc -l)
echo "Found files to process: $total_files + gp: $gp_files"
#------------------------------------------------------

counter=0

echo "Processing non-METIS files from $routine_dir"
# Iterate through each .out file in the method directory
find "$routine_dir" -name '*.out' | grep "${routine}_" | while read -r file_path; do
  process_file $file_path
  counter=$(($counter + 1))
  show_progress $counter $total_files
done

echo "Processing METIS files from $routine_gp_dir"
find "$routine_gp_dir" -name '*.out' | grep "${routine}" | while read -r file_path; do
  process_file $file_path
  counter=$(($counter + 1))
  show_progress $counter $gp_files
done

echo "Processed $counter files"
echo "Results saved to $output_dir"
