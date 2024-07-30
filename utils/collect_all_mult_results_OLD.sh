#!/bin/bash

#************************************
#*************  INPUT  **************
#************************************************************************************************************************************************

# Default values
root_dir="../../../Downloads/outputs_2024-07-26/SbatchMan/outputs/marzola"
routine="spmmcsr"
method="ALL" # Options: saad, metis, denseAMP, original, clubs, ALL

# Function to display usage
usage() {
    echo "Usage: $0 [-r root_dir] [-u routine] [-m method]"
    echo "  -r  Root directory (default: ../../../Downloads/mult_outputs_24)"
    echo "  -u  Routine (default: spmvcsr)"
    echo "  -m  Method (default: ALL)"
    exit 1
}

# Parse options
while getopts ":r:u:m:" opt; do
    case ${opt} in
        r )
            root_dir=$OPTARG
            ;;
        u )
            routine=$OPTARG
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

# Print the values to check
echo "Root directory: $root_dir"
echo "Routine: $routine"
echo "Routine directory: $routine_dir"
echo "Routine GP directory: $routine_gp_dir"
echo "Method: $method"

if [ ! -d "$routine_dir" ]; then
  echo "ERROR: COULD NOT FIND RESULT DIRECTORY $routine_dir"
  exit 1
fi

# Output CSV file
output_dir="results/results_2024/mult/$routine"
mkdir -p "$output_dir"

#************************************
#************************************
#************************************************************************************************************************************************

#************************************************************************************************************************************************
#************* UTILS    *************
#************************************

get_outfile() {
  local algo=$1
  echo "$output_dir/$algo/results_${routine}_$algo.csv"
}

get_header() {
  local algo=$1
  case $algo in
    clubs)
      echo "routine matrix algo mask centroid tau rows cols nnz time"
      ;;
    metis-edge-cut)
      echo "routine matrix algo objective parts rows cols nnz time"
      ;;
    metis-volume)
      echo "routine matrix algo objective parts rows cols nnz time"
      ;;
    saad)
      echo "routine matrix algo tau rows cols nnz time"
      ;;
    denseAMP)
      echo "routine matrix algo tau rows cols nnz time"
      ;;
    original)
      echo "routine matrix algo rows cols nnz time"
      ;;
    *)
      exit 1
      ;;
  esac
}

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

# Function to extract time from a file
extract_time_and_variables() {
  local file_path=$1

  # Extract time
  local time=$(grep 'calculation finished in' "$file_path" | head -1 | awk '{print $(NF-1)}')

  # Extract mb, nb, nnzb
  if [[ $routine == *"bsr"* ]]; then
    m=$(grep '^mb =' "$file_path" | awk '{print $3}')
    n=$(grep '^nb =' "$file_path" | awk '{print $3}')
    nnz=$(grep '^nnzb =' "$file_path" | awk '{print $3}')
  else
    m=$(grep '^m =' "$file_path" | awk '{print $3}')
    n=$(grep '^n =' "$file_path" | awk '{print $3}')
    nnz=$(grep '^nnz =' "$file_path" | awk '{print $3}')
  fi
  echo "$time" "$m" "$n" "$nnz"
}
# Function to extract variables from clubs output
extract_variables_club() {
  local filename=$1
  IFS='_' read -ra PARTS <<< "$filename"
  local matrix=$( echo ${PARTS[1]} | sed "s/-mtx"// )
  local collection=${PARTS[2]}
  collection=${collection#"$matrix-"}
  IFS='-' read -ra VARS <<< "$collection"

  local cents=$(echo ${VARS[0]} | sed 's/cent//')
  local msk=$(echo ${VARS[1]} | sed 's/msk//')
  local tau=$(echo ${VARS[2]} | sed 's/tau//')

  echo "$matrix clubs $msk $cents $tau"
}

# Function to extract variables from saad output
extract_variables_saad() {
  local filename=$1
  IFS='_' read -ra PARTS <<< "$filename"
  local matrix=$( echo ${PARTS[1]} | sed "s/-mtx"// )
  local collection=${PARTS[2]}
  collection=${collection#"$matrix-"}
  IFS='-' read -ra VARS <<< "$collection"

  local tau=$(echo ${VARS[1]} | sed "s/tau"// )

  echo "$matrix saad $tau"
}

extract_variables_denseAMP() {
  local filename=$1
  IFS='_' read -ra PARTS <<< "$filename"
  local matrix=$( echo ${PARTS[1]} | sed "s/-mtx"// )
  local collection=${PARTS[2]}
  collection=${collection#"$matrix-"}
  IFS='-' read -ra VARS <<< "$collection"

  local tau=$(echo ${VARS[1]} | sed "s/tau"// )

  echo "$matrix denseAMP $tau"
}

extract_variables_metis() {
  local filename=$1
  IFS='_' read -ra PARTS <<< "$filename"
  local matrix=$( echo ${PARTS[1]} | sed "s/-mtx"// )
  local collection=${PARTS[2]}
  collection=${collection#"$matrix-metis"}
  local obj
  local parts
  IFS='-' read -ra VARS <<< "$collection"

  if [[ $collection == *"volume"* ]]; then
    obj="volume"
    parts=$(echo ${VARS[2]} | sed 's/part//')
  else
    obj="edge-cut"
    parts=$(echo ${VARS[3]} | sed 's/part//')
  fi

  #TODO ADD PARAMS
  echo "$matrix metis $obj $parts"
}

# Function to extract variables from original output
extract_variables_original() {
  local filename=$1
  IFS='_' read -ra PARTS <<< "$filename"
  local matrix=$( echo ${PARTS[1]} | sed "s/-mtx"// )

  echo "$matrix original"
}


find_file_type() {
  local file_name=$1

  if [[ $file_name == *"msk"* && $file_name == *"cent"* && $file_name == *"tau"* ]]; then
    file_type="clubs"
  elif [[ $file_name == *"metis"* ]]; then # TODO CHECK
    if [[ $file_name == *"edge"* ]];then
      file_type="metis-edge-cut"
    else
      file_type="metis-volume"
    fi
  elif [[ $file_name == *"Saad"* ]]; then # TODO CHECK
    file_type="saad"
  elif [[ $file_name == *"DenseAMP"* ]]; then # TODO CHECK
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

  local time_vars=$(extract_time_and_variables "$file_path")

  local time=$(echo "$time_vars" | awk '{print $1}')
  local m=$(echo "$time_vars" | awk '{print $2}')
  local n=$(echo "$time_vars" | awk '{print $3}')
  local nnz=$(echo "$time_vars" | awk '{print $4}')

  file_type=$(find_file_type "$file_name")
  if [[ "$method" != "$file_type" && "$method" != "ALL" ]]; then
     return 0
  fi

  case $file_type in
    clubs)
      variables=$(extract_variables_club "$file_name")
      ;;
    metis-edge-cut)
      variables=$(extract_variables_metis "$file_name")
      ;;
    metis-volume)
      variables=$(extract_variables_metis "$file_name")
      ;;
    saad)
      variables=$(extract_variables_saad "$file_name")
      ;;
    denseAMP)
      variables=$(extract_variables_denseAMP "$file_name")
      ;;
    original)
      variables=$(extract_variables_original "$file_name")
      ;;
    *)
      echo "COULD NOT CLASSIFY $file_name"
      return 1
      ;;
  esac

  #echo "--Assigned to $file_type"
  #echo "--printing results in $(get_outfile "$file_type")"
  echo "$routine $variables $m $n $nnz $time" >> "$(get_outfile "$file_type")"
}


#************************************
#************************************
#************************************************************************************************************************************************

#************************************************************************************************************************************************

#creates out files for processed methods
if [[ "$method" == "ALL" ]]; then
  algos=( "clubs" "metis-edge-cut" "metis-volume" "saad" "original" "denseAMP")
elif [[ "$method" == "metis" ]]; then
  algos=( "metis-edge-cut" "metis-volume")
else
  algos=( "$method" )
fi

for algo in "${algos[@]}"; do
  mkdir -p "$output_dir/$algo"
  echo "creating csv:" "$(get_outfile "$algo")"
  get_header "$algo" > "$(get_outfile "$algo")" # Print header for each algo
done
#-----------------------------------------



#find files
total_files=$(find "$routine_dir" -name '*.out' | wc -l)
gp_files=$(find "$routine_gp_dir" -name '*.out' | wc -l)
echo "Found files to process: $total_files + gp: $gp_files"
#------------------------------------------------------

counter=0

if [[ "$method" != "metis-edge-cut" && "$method" != "metis-volume" ]]; then
    echo "Processing non-METIS files"
  # Iterate through each .out file in the method directory
  find "$routine_dir" -name '*.out' | grep "${routine}_" | while read -r file_path; do
    process_file $file_path
    counter=$(($counter + 1))
    show_progress $counter $total_files
  done
fi

if [[ "$method" == "metis-edge-cut" || "$method" == "metis-volume" || "$method" == "ALL" ]]; then
  echo "Processing METIS files"
  find "$routine_gp_dir" -name '*.out' | grep "${routine}" | while read -r file_path; do
    process_file $file_path
    counter=$(($counter + 1))
    show_progress $counter $gp_files
  done
fi

echo "Processed $counter files"
echo "Results saved to $output_dir"
