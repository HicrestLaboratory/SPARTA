#!/bin/bash

#************************************
#*************  INPUT  **************
#************************************************************************************************************************************************

# Default values
root_dir="results/results_2024/mult_data"
csv_dir="results/results_2024/mult_csv"
routine="spmmcsr"
method="ALL" # Options: saad, metis, denseAMP, original, clubs, ALL

# Function to display usage
usage() {
    echo "Usage: $0 [-r root_dir] [-u routine] [-m method]"
    echo "  -r  Root directory (default: $root_dir)"
    echo "  -u  Routine (default: $routine)"
    echo "  -m  Method (default: $method)"
    echo "  -c  CSV directory (default $csv_dir)"
    exit 1
}

# Parse options
while getopts ":r:u:m:c:" opt; do
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
        c )
            csv_dir=$OPTARG
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
echo "root directory: $root_dir"
echo "Routine: $routine"
echo "Method: $method"

if [ ! -d "$root_dir" ]; then
  echo "ERROR: COULD NOT FIND RESULT DIRECTORY $root_dir"
  exit 1
fi

# Output CSV file
output_dir="$csv_dir/$routine"
mkdir -p "$output_dir"
echo "storing csvs in output directory: $output_dir"

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
  local matrix=${PARTS[1]%-mtx}
  local collection=${PARTS[2]#"$matrix-"}
  IFS='-' read -ra VARS <<< "$collection"

  local cents=${VARS[0]#cent}
  local msk=${VARS[1]#msk}
  local tau=${VARS[2]#tau}

  echo "$matrix clubs $msk $cents $tau"
}

# Function to extract variables from saad output
extract_variables_saad() {
  local filename=$1
  IFS='_' read -ra PARTS <<< "$filename"
  local matrix=${PARTS[1]%-mtx}
  local collection=${PARTS[2]#"$matrix-"}
  IFS='-' read -ra VARS <<< "$collection"

  local tau=${VARS[1]#tau}

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
  local matrix=${PARTS[1]%-mtx}
  local collection=${PARTS[2]#"$matrix-metis"}
  local obj
  local parts
  IFS='-' read -ra VARS <<< "$collection"

  if [[ $collection == *"volume"* ]]; then
    obj="volume"
    parts=${VARS[2]#part}
  else
    obj="edge-cut"
    parts=${VARS[3]#part}
  fi

  echo "$matrix metis $obj $parts"
}

# Function to extract variables from original output
extract_variables_original() {
  local filename=$1
  IFS='_' read -ra PARTS <<< "$filename"
  local matrix=$( echo ${PARTS[1]} | sed "s/-mtx"// )

  echo "$matrix original"
}


process_file()
{
  local file_path=$1
  local algo=$2

  local file_name=$(basename "$file_path")

  local time_vars=$(extract_time_and_variables "$file_path")

  local time=$(echo "$time_vars" | awk '{print $1}')
  local m=$(echo "$time_vars" | awk '{print $2}')
  local n=$(echo "$time_vars" | awk '{print $3}')
  local nnz=$(echo "$time_vars" | awk '{print $4}')

  case $algo in
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
      echo "COULD NOT CLASSIFY $file_name for $algo"
      return 1
      ;;
  esac

  #echo "--Assigned to $file_type"
  #echo "--printing results in $(get_outfile "$file_type")"
  echo "$routine $variables $m $n $nnz $time" >> "$(get_outfile "$algo")"
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

for algo in "${algos[@]}";do
    counter=0
    algo_dir=$root_dir/$routine/$algo
    total_files=$(find "$algo_dir" -name '*.out' | wc -l)
    echo "PROCESSING $total_files mult files from $algo_dir for $algo"
    find "$algo_dir" -name '*.out' | grep "${routine}" | while read -r "file_path"; do
      process_file $file_path $algo
      counter=$(($counter + 1))
      if (( counter % 10 == 0 )); then
        show_progress $counter $total_files
      fi
    done
    echo "Processed $counter files"
done

echo "Results saved to $output_dir"
