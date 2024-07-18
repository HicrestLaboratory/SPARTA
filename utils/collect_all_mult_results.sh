#!/bin/bash

# Define the root directory
root_dir="Downloads/club_outputs/SbatchMan/outputs/marzola"

# Output CSV file
output_file="extracted_mult_times.csv"
method="spmmcsr"
method_dir="$root_dir/$method"


# Initialize the CSV file with headers
echo "method,structure,mask,centroid,tau,time" > $output_file

# Function to extract time from a file
extract_time() {
  local file_path=$1
  grep 'calculation finished in' $file_path | head -1 | awk '{print $5}'
}

# Function to extract variables from clubs outout
extract_variables_structure1() {
  local filename=$1
  IFS='-' read -ra PARTS <<< "$filename"
  local mask=$(echo ${PARTS[2]} | sed 's/msk//')
  local centroid=$(echo ${PARTS[1]} | sed 's/cent//')
  local tau=$(echo ${PARTS[3]} | sed 's/tau//')
  echo "clubs,$msk,$cent,$tau"
}

# Function to extract variables from original output TODO
extract_variables_structure2() {
  local filename=$1
  IFS='-' read -ra PARTS <<< "$filename"
  echo "original"
}

if [ -d "$method_dir" ]; then
    # Iterate through each .out file in the method directory
    find "$method_dir" -name '*.out' -o -name '*.err' | while read file_path; do
    time=$(extract_time $file_path)
    if [ ! -z "$time" ]; then
        file_name=$(basename $file_path)
        if [[ $file_name == *"msk"* && $file_name == *"cent"* && $file_name == *"tau"* ]]; then
        variables=$(extract_variables_structure1 $file_name)
        else
        variables=$(extract_variables_structure2 $file_name)
        fi
        echo "$method,$variables,$time" >> $output_file
    fi
    done
fi

echo "Results saved to $output_file"
