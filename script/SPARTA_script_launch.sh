#!/bin/bash

#SBATCH --job-name=sparta-test
#SBATCH --output=sparta-test-%j.out
#SBATCH --error=sparta-test-%j.err
#SBATCH --time=100:00:00
#SBATCH --partition training
#SBATCH --gres=gpu
#SBATCH --nodes 1
#SBATCH --ntasks=1


module load cuda-11.4.0
cd ${HOME}/SPARTA
params="";
for i in "$@"; 
   do
	params=" ${params} $d $i"
   done
${params} #call any script from here
