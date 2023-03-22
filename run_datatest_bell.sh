#!/bin/bash

path=../SPARTA_datasets/suitsparse_collection_3

Matrices=('ckt11752_dc_1.mtx' 'email-Enron.mtx' 'FEM_3D_thermal1.mtx' 'g7jac100sc.mtx' 'k3plates.mtx' 'lhr17.mtx' 'poisson3Da.mtx' 'TEM27623.mtx' 'vsp_south31_slptsk.mtx')

for i in ${!Matrices[@]}
do
    for script in batch/BELLPACK
    do
        for blk_size in 32 64 128 256 512 1024
        do
            for Bcols in 2048 8192
            do
                echo "sbatch $script $path/${Matrices[$i]} $blk_size $Bcols"
                sbatch $script $path/${Matrices[$i]} $blk_size $Bcols
            done
        done
    done
done
