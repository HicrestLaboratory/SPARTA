#!/bin/bash

#SBATCH --job-name=sparta-test
#SBATCH --output=sparta-test-%j.out
#SBATCH --error=sparta-test-%j.err

#SBATCH --partition training
#SBATCH --gres=gpu
#SBATCH --nodes 1
#SBATCH --time=00:01:00

#SBATCH --ntasks=1


export OPTS="-w 1 -r 5 -v -1 -i 4"
module load cuda-10.2
module load gcc-9.3.0
cd /home/clusterusers/pasyloslabini/SPARTA
make clean
make test_cublas_cusparse_comparison
cd /home/clusterusers/pasyloslabini/SPARTA
source /home/clusterusers/pasyloslabini/SPARTA/script/cusparse-vs-gemm-test.sh
