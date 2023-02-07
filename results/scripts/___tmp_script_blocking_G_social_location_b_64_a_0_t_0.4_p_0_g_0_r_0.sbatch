#!/bin/bash -l
#SBATCH --job-name=blocking_G_social_location_b_64_a_0_t_0.4_p_0_g_0_r_0
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --output=results//scripts/outputs/blocking_G_social_location_b_64_a_0_t_0.4_p_0_g_0_r_0.%j.o
#SBATCH --error=results//scripts/errors/blocking_G_social_location_b_64_a_0_t_0.4_p_0_g_0_r_0.%j.e
#SBATCH --account=g34
#SBATCH --partition=normal
#SBATCH --constraint=ssd

./programs/general/Matrix_Blocking -f data/real_world//social_location.el -b 64 -t 0.4 -a 0 -r 0 -p 0 -g 0 -v 1 -o results//social_location/blocking_G_social_location_b_64_a_0_t_0.4_p_0_g_0_r_0.txt -P 1 -n blocking_G_social_location_b_64_a_0_t_0.4_p_0_g_0_r_0

sleep 1s

