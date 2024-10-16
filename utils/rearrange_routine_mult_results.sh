#!/bin/bash

# Default values
root_dir="../../../Downloads/outputs_02_10/SbatchMan/outputs/marzola"
routine="spmmcsr"
collected_data_folder="results/results_$(date +'%d_%m_%Y')/mult_data"
clean_folders="0"
method="all"

# Function to display usage
usage() {
    echo "Usage: $0 [-r root_dir] [-u routine] [-x clean folders]"
    echo "  -r  Root directory (default: $root_dir)"
    echo "  -o  Data storage directory (default: $collected_data_folder)"
    echo "  -u  Routine (default: $routine)"
    echo "  -x  Clean output folders"
    echo "  -m Method (metis-edge-cut, metis-volume, clubs, denseAMP, saad, original, patoh)"
    exit 1
}

# Parse options
while getopts ":r:u:o:m:x" opt; do
    case ${opt} in
        r )
            root_dir=$OPTARG
            ;;
        u )
            routine=$OPTARG
            ;;
        x )
            clean_folders="1"
            ;;
        o )
            collected_data_folder=$OPTARG
            ;;
        m )
            method=$OPTARG
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
routine_patoh_dir="$root_dir/${routine}patoh/finished"

# Print the values to check
echo "Root directory: $root_dir"
echo "Data collection folder: $collected_data_folder"
echo "Routine: $routine"
echo "Routine directory: $routine_dir"
echo "Routine GP directory: $routine_gp_dir"

if [ ! -d "$routine_dir" ]; then
  echo "ERROR: COULD NOT FIND RESULT DIRECTORY $routine_dir"
  exit 1
fi

out_routine_dir=$collected_data_folder/$routine
echo $out_routine_dir
mkdir -p "$out_routine_dir"

if [ "$method" != "all" ]; then
  algos=( "$method" )
else
  algos=( "metis-edge-cut" "metis-volume" "clubs" "denseAMP" "saad" "original" "patoh" )
fi

for algo in "${algos[@]}"; do
  if [ "$clean_folders" -eq 1 ]; then
    rm -rf "$out_routine_dir/$algo"
  fi
  mkdir -p "$out_routine_dir/$algo"
done

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
  elif [[ $file_name == *"PaToH"* ]]; then
    file_type="patoh"
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

# Find files
total_files=$(find "$routine_dir" -name '*.out' | wc -l)
gp_files=$(find "$routine_gp_dir" -name '*.out' | wc -l)
patoh_files=$(find "$routine_patoh_dir" -name '*.out' | wc -l)
echo "Found files to process: $total_files + gp: $gp_files + patoh $patoh_files"

counter=0

# Process non-METIS files
if [[ " ${algos[@]} " =~ " clubs " ]]; then
  echo "Processing non-METIS files from $routine_dir"
  find "$routine_dir" -name '*.out' | grep "${routine}_" | while read -r file_path; do
    process_file "$file_path"
    counter=$((counter + 1))
    show_progress $counter $total_files
  done
fi

# Process METIS files
if [[ " ${algos[@]} " =~ "metis-edge-cut" ]]; then
  echo "Processing METIS files from $routine_gp_dir"
  find "$routine_gp_dir" -name '*.out' | grep "${routine}" | while read -r file_path; do
    process_file "$file_path"
    counter=$((counter + 1))
    show_progress $counter $gp_files
  done
fi

# Process PATOH files
if [[ " ${algos[@]} " =~ " patoh " ]]; then
  echo "Processing PATOH files from $routine_patoh_dir"
  find "$routine_patoh_dir" -name '*.out' | grep "${routine}" | while read -r file_path; do
    process_file $file_path
    counter=$(($counter + 1))
    show_progress $counter $patoh_files
  done
fi