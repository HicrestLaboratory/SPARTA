#!/bin/bash -l
#SBATCH --job-name=blocking_G_cs_department_b_16_a_2_B_32_r_0
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --output=results//scripts/outputs/blocking_G_cs_department_b_16_a_2_B_32_r_0.%j.o
#SBATCH --error=results//scripts/errors/blocking_G_cs_department_b_16_a_2_B_32_r_0.%j.e
#SBATCH --account=g34
#SBATCH --partition=normal
#SBATCH --constraint=ssd

./programs/general/Matrix_Blocking -f data/real_world//cs_department.el -b 16 -a 2 -B 32 -r 0 -v 1 -o results//cs_department/blocking_G_cs_department_b_16_a_2_B_32_r_0.txt -P 1 -n blocking_G_cs_department_b_16_a_2_B_32_r_0

sleep 1s

