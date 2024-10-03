#!/bin/bash

initout() {
        current_date=$(date +%Y-%m-%d)  # Format: YYYY-MM-DD
        current_time=$(date +%H:%M:%S)

        echo "# ----- Init file ${current_date} ${current_time} -----"
}

runcommand() {
        commandToRun=$1
        current_date=$(date +%Y-%m-%d)  # Format: YYYY-MM-DD
        current_time=$(date +%H:%M:%S)

        echo "# ====================================================="
        echo "Date: ${current_date} Time: ${current_time}"
        echo "Running ${commandToRun}..."
        echo "# -----------------------------------------------------"
        ${commandToRun}
        current_date=$(date +%Y-%m-%d)  # Format: YYYY-MM-DD
        current_time=$(date +%H:%M:%S)
        echo "# -----------------------------------------------------"
        echo "Done"
        echo "Date: ${current_date} Time: ${current_time}"
        echo "# ====================================================="
}

my_dataset="datasetA"
#my_dataset="toy"

current_date=$(date +%Y-%m-%d)  # Format: YYYY-MM-DD
current_time=$(date +%H:%M:%S)

initout
runcommand "./utils/collect_all_reord_results.sh -f /data/clubs/Club_dataset/${my_dataset} -r /data/clubs/Club_result/${my_dataset} -o results_${current_date} -s patoh" 
runcommand "tar -cf results_${current_date}.tar results_${current_date}/"
runcommand "gzip results_${current_date}.tar"
