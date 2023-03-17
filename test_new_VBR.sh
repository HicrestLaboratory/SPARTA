#!/bin/bash
for Matrix in data/minitest/*.el
do
    for script in batch/M5a5F1
    do
        for blk_size in 128
        do
            for Bcols in 8192
            do
                echo "sbatch $script $Matrix $blk_size $Bcols"
                sbatch $script $Matrix $blk_size $Bcols 8
            done
        done
    done
done
