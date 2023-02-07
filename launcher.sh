#!/bin/bash -l
#SBATCH --job-name="reorder_psylos"
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --output=sparta.%j.o
#SBATCH --error=sparta.%j.e
#SBATCH --account="g34"
#SBATCH --partition=normal
#SBATCH --constraint=ssd

./src/scripts/run_blocking_experiments.sh data/real_world/ results/ programs/general/Matrix_Blocking
