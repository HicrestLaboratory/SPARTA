#!/bin/bash

#SBATCH --job-name=sparta-test
#SBATCH --output=sparta-test-%j.out
#SBATCH --error=sparta-test-%j.err
#SBATCH --time=04:00:00
#SBATCH --partition dgx2
#SBATCH -n 1

#SBATCH -G 1


cd ${HOME}/SPARTA
#ls 
#nvidia-smi
source ${1}  #call any script from here
