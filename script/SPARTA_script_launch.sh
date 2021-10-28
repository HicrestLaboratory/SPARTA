#!/bin/bash

#SBATCH --job-name=sparta-test
#SBATCH --output=sparta-test-%j.out
#SBATCH --error=sparta-test-%j.err
#SBATCH --time=04:00:00
#SBATCH --partition amda100
#SBATCH --nodes 1

#SBATCH --ntasks=1


module load cuda
cd ${HOME}/SPARTA
#ls 
#nvidia-smi
source ${1}  #call any script from here
