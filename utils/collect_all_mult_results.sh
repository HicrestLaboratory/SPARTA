#!/bin/bash

# Define the root directory
root_dir="../../../Downloads/mult_outputs"
routine="spmvbsr"
routine_dir="$root_dir/$routine/finished"

if [ ! -d "$routine_dir" ]; then
  echo "ERROR: COULD NOT FIND RESULT DIRECTORY $routine_dir"
  exit 1
fi

# Output CSV file
output_dir="results/results_2024/mult/$routine"
mkdir -p "$output_dir"

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
    metis)
      echo "routine matrix algo objective parts rows cols nnz time"
      ;;
    saad)
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

  #TODO ADD TAU
  echo "$matrix saad"
}

extract_variables_metis() {
  local filename=$1
  IFS='_' read -ra PARTS <<< "$filename"
  local matrix=$( echo ${PARTS[1]} | sed "s/-mtx"// )
  local collection=${PARTS[2]}
  collection=${collection#"$matrix-"}
  IFS='-' read -ra VARS <<< "$collection"

  #TODO ADD PARAMS

  collection=${collection#"$matrix-"}
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
    file_type="metis"
  elif [[ $file_name == *"saad"* ]]; then # TODO CHECK
    file_type="saad"
  else
    file_type="original"
  fi

  echo $file_type
}

process_file()
{
  local file_path=$1

  local file_name=$(basename "$file_path")
  echo "Processing: $file_name"

  local time_vars=$(extract_time_and_variables "$file_path")

  local time=$(echo "$time_vars" | awk '{print $1}')
  local m=$(echo "$time_vars" | awk '{print $2}')
  local n=$(echo "$time_vars" | awk '{print $3}')
  local nnz=$(echo "$time_vars" | awk '{print $4}')

  file_type=$(find_file_type "$file_name")
  case $file_type in
    clubs)
      echo "$(extract_variables_club "$file_name")"
      variables=$(extract_variables_club "$file_name")
      ;;
    metis)
      variables=$(extract_variables_metis "$file_name")
      ;;
    saad)
      variables=$(extract_variables_saad "$file_name")
      ;;
    original)
      variables=$(extract_variables_original "$file_name")
      ;;
    *)
      exit 1
      ;;
  esac

  #echo "--Assigned to $file_type"
  #echo "--printing results in $(get_outfile "$file_type")"
  echo "$routine $variables $m $n $nnz $time" >> "$(get_outfile "$file_type")"
}

algos=( "clubs" "metis" "saad" "original" )
for algo in "${algos[@]}"; do
  mkdir -p "$output_dir/$algo"
  echo "creating csv:" "$(get_outfile "$algo")"
  get_header "$algo" > "$(get_outfile "$algo")" # Print header for each algo
done


total_files=$(find "$routine_dir" -name '*.out' | wc -l)
echo "Found files to process: $total_files"


counter=0
# Iterate through each .out file in the method directory
find "$routine_dir" -name '*.out' | grep "${routine}_" | while read -r file_path; do
  process_file $file_path
  counter=$((counter + 1))
done

echo "Processed $counter files out of $total_files"
echo "Results saved to $output_dir"
