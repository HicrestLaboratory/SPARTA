#!/bin/bash

sbatch batch/exp_1 data/real_world/wiki-Vote_r.el 16
sbatch batch/exp_1 data/real_world/wiki-Vote_r.el 32
sbatch batch/exp_1 data/real_world/wiki-Vote_r.el 64

sbatch batch/exp_2 data/real_world/wiki-Vote_r.el 16
sbatch batch/exp_2 data/real_world/wiki-Vote_r.el 32
sbatch batch/exp_2 data/real_world/wiki-Vote_r.el 64

sbatch batch/exp_3 data/real_world/wiki-Vote_r.el 16
sbatch batch/exp_3 data/real_world/wiki-Vote_r.el 32
sbatch batch/exp_3 data/real_world/wiki-Vote_r.el 64

sbatch batch/exp_5 data/real_world/wiki-Vote_r.el 16
sbatch batch/exp_5 data/real_world/wiki-Vote_r.el 32
sbatch batch/exp_5 data/real_world/wiki-Vote_r.el 64

# --------------------------------------------------

sbatch batch/exp_1 data/real_world/wiki-Vote_r.el 16 2048
sbatch batch/exp_1 data/real_world/wiki-Vote_r.el 32 2048
sbatch batch/exp_1 data/real_world/wiki-Vote_r.el 64 2048

sbatch batch/exp_2 data/real_world/wiki-Vote_r.el 16 2048
sbatch batch/exp_2 data/real_world/wiki-Vote_r.el 32 2048
sbatch batch/exp_2 data/real_world/wiki-Vote_r.el 64 2048

sbatch batch/exp_3 data/real_world/wiki-Vote_r.el 16 2048
sbatch batch/exp_3 data/real_world/wiki-Vote_r.el 32 2048
sbatch batch/exp_3 data/real_world/wiki-Vote_r.el 64 2048

sbatch batch/exp_5 data/real_world/wiki-Vote_r.el 16 2048
sbatch batch/exp_5 data/real_world/wiki-Vote_r.el 32 2048
sbatch batch/exp_5 data/real_world/wiki-Vote_r.el 64 2048


# --------------------------------------------------

sbatch batch/exp_1 data/real_world/wiki-Vote_r.el 16 4096
sbatch batch/exp_1 data/real_world/wiki-Vote_r.el 32 4096
sbatch batch/exp_1 data/real_world/wiki-Vote_r.el 64 4096

sbatch batch/exp_2 data/real_world/wiki-Vote_r.el 16 4096
sbatch batch/exp_2 data/real_world/wiki-Vote_r.el 32 4096
sbatch batch/exp_2 data/real_world/wiki-Vote_r.el 64 4096

sbatch batch/exp_3 data/real_world/wiki-Vote_r.el 16 4096
sbatch batch/exp_3 data/real_world/wiki-Vote_r.el 32 4096
sbatch batch/exp_3 data/real_world/wiki-Vote_r.el 64 4096

sbatch batch/exp_5 data/real_world/wiki-Vote_r.el 16 4096
sbatch batch/exp_5 data/real_world/wiki-Vote_r.el 32 4096
sbatch batch/exp_5 data/real_world/wiki-Vote_r.el 64 4096
