#!/bin/bash -l
#SBATCH --job-name=blocking_G_ia-wikiquote-user-edits-nodup_b_256_a_1_t_0.7_p_1_g_0_r_0
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --output=results//scripts/outputs/blocking_G_ia-wikiquote-user-edits-nodup_b_256_a_1_t_0.7_p_1_g_0_r_0.%j.o
#SBATCH --error=results//scripts/errors/blocking_G_ia-wikiquote-user-edits-nodup_b_256_a_1_t_0.7_p_1_g_0_r_0.%j.e
#SBATCH --account=g34
#SBATCH --partition=normal
#SBATCH --constraint=ssd

./programs/general/Matrix_Blocking -f data/real_world//ia-wikiquote-user-edits-nodup.el -b 256 -t 0.7 -a 1 -r 0 -p 1 -g 0 -v 1 -o results//ia-wikiquote-user-edits-nodup/blocking_G_ia-wikiquote-user-edits-nodup_b_256_a_1_t_0.7_p_1_g_0_r_0.txt -P 1 -n blocking_G_ia-wikiquote-user-edits-nodup_b_256_a_1_t_0.7_p_1_g_0_r_0

sleep 1s

