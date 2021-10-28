#!/bin/bash

#SBATCH --job-name=sparta-test
#SBATCH --output=sparta-test-%j.out
#SBATCH --error=sparta-test-%j.err

#SBATCH --partition training
#SBATCH --gres=gpu
#SBATCH --nodes 1
#SBATCH --time=99:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --ntasks=1


export OPTS="-w 1 -r 5 -v 1 -i 4"
module load cuda
cd /home/clusterusers/pasyloslabini/SPARTA
make clean
make test_cublas_VBS
cd /home/clusterusers/pasyloslabini/SPARTA
./programs/cuda/test_cublas_VBS ${OPTS} -m 320000 -n 128 -k 320000 -p 128 -b 0.01 -q 0.06
