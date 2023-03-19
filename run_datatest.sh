#!/bin/bash

if [[ $# < 1 ]]
then
    echo "usage: ./run_dataset path_to_dataset"
else

dataset=$1

for Matrix in $dataset/*.mtx
do
    for script in batch/VBR batch/BELLPACK
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

    for script in batch/CSR
    do
        for Bcols in 1024 2048 4096 8192
        do
            echo "sbatch $script $Matrix $Bcols"
            sbatch $script $Matrix $Bcols
        done
    done

done

fi
