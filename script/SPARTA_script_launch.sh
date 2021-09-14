#!/bin/bash

#SBATCH --job-name=sparta-test
#SBATCH --output=sparta-test-%j.out
#SBATCH --error=sparta-test-%j.err

#SBATCH --partition training
#SBATCH --gres=gpu
#SBATCH --nodes 1
#SBATCH --time=99:00:00

#SBATCH --ntasks=1


module load cuda-11.4.0
module load gcc-9.3.0
cd /home/clusterusers/pasyloslabini/SPARTA
source ${1}  #call any script from here
