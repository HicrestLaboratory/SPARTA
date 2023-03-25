#!/bin/bash

path=../SPARTA_datasets/suitsparse_collection_3

Matrices=('ckt11752_dc_1.mtx' 'email-Enron.mtx' 'FEM_3D_thermal1.mtx' 'g7jac100sc.mtx' 'k3plates.mtx' 'lhr17.mtx' 'poisson3Da.mtx' 'TEM27623.mtx' 'vsp_south31_slptsk.mtx')
Taus=(0.2  0.2  0.2  0.2 0.2  0.2  0.2  0.2  0.2)


for blk_size in 32 64 128 256 512 1024
do
    echo "blk_size = $blk_size"

    for i in ${!Matrices[@]}
    do
#       echo "With $row_blk_size x $col_blk_size, matrix ${Matrices[$i]} has tau ${Taus[$i]}"

      for c in 2048 8192
      do

        if [ -f "$path/${Matrices[$i]}" ]
        then

#         echo "sbatch batch/BELLPACK $path/${Matrices[$i]} $blk_size ${Taus[$i]}"
        sbatch batch/BELLPACK $path/${Matrices[$i]} $blk_size ${Taus[$i]} $c

#         echo "sbatch batch/CUTLASS_BELLPACK $path/${Matrices[$i]} $blk_size ${Taus[$i]}"
        sbatch batch/CUTLASS_BELLPACK $path/${Matrices[$i]} $blk_size ${Taus[$i]} $c

        else
          echo "ERROR: $path/${Matrices[$i]} does not exists."
        fi

      done
    done

done
