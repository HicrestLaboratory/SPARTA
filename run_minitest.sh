#!/bin/bash
for Matrix in data/minitest/*.el
do
    for script in batch/csrMult batch/BellFix batch/BellMult batch/VbrFix batch/VbrNoFix
    do
        for blk_size in 16 32 64 128
        do
            for Bcols in 1024 2048 4096 8192
            do
                echo "sbatch $script $Matrix $blk_size $Bcols"
                sbatch $script $Matrix $blk_size $Bcols
            done
        done
    done
done
