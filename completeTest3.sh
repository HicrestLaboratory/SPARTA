#!/bin/bash

path=../SPARTA_datasets/suitsparse_big

Matrices=('2cubes_sphere.mtx' 'boneS01.mtx' 'engine.mtx' 'G_n_pin_pout.mtx' 'loc-Gowalla.mtx' 'n4c6-b9.mtx' 'rgg_n_2_17_s0.mtx' 'vsp_finan512_scagr7-2c_rlfddd.mtx' '598a.mtx' 'c-73.mtx' 'FEM_3D_thermal2.mtx' 'Goodwin_095.mtx' 'matrix_9.mtx' 'para-10.mtx' 'Si41Ge41H72.mtx' 'barrier2-9.mtx' 'crashbasis.mtx' 'fullb.mtx' 'Goodwin_127.mtx' 'mono_500Hz.mtx' 'para-7.mtx' 'torso1.mtx' 'bmwcra_1.mtx' 'Dubcova3.mtx' 'Ga10As10H30.mtx' 'hvdc2.mtx' 'n4c6-b7.mtx' 'pkustk14.mtx' 'twotone.mtx')
Taus=( 0.1  0.1  0.1  0.1  0.1 0.1  0.1  0.1  0.1 0.1  0.1  0.1  0.1  0.1 0.1  0.1  0.1  0.1 0.1 0.1  0.1  0.1  0.1  0.1 0.1  0.1  0.1  0.1 0.1 )

for row_blk_size in 32 64 128 256 512 1024
do
    col_blk_size=$row_blk_size

    for i in ${!Matrices[@]}
    do
        echo "With $row_blk_size x $row_blk_size, matrix ${Matrices[$i]} has tau ${Taus[$i]}"

        if [ -f "$path/${Matrices[$i]}" ]
        then

#           echo "sbatch batch/VBR_rect_a5 $path/${Matrices[$i]} $row_blk_size $row_blk_size ${Taus[$i]}"
          sbatch batch/VBR_rect_a5 $path/${Matrices[$i]} $row_blk_size $row_blk_size ${Taus[$i]}

#           echo "sbatch batch/VBR_rect_a2 $path/${Matrices[$i]} $row_blk_size $row_blk_size"
          sbatch batch/VBR_rect_a2 $path/${Matrices[$i]} $row_blk_size $row_blk_size

          #           echo "sbatch batch/VBR_rect_a5 $path/${Matrices[$i]} $row_blk_size $row_blk_size ${Taus[$i]}"
          sbatch batch/VBR_batch_a5 $path/${Matrices[$i]} $row_blk_size $row_blk_size ${Taus[$i]}

#           echo "sbatch batch/VBR_rect_a2 $path/${Matrices[$i]} $row_blk_size $row_blk_size"
          sbatch batch/VBR_batch_a2 $path/${Matrices[$i]} $row_blk_size $row_blk_size

          if [ $row_blk_size == $row_blk_size ]
          then
#           echo "sbatch batch/BELLPACK_a5 $path/${Matrices[$i]} $row_blk_size ${Taus[$i]}"
          sbatch batch/BELLPACK_a5 $path/${Matrices[$i]} $row_blk_size ${Taus[$i]}

          #           echo "sbatch batch/BELLPACK_a2 $path/${Matrices[$i]} $row_blk_size"
          sbatch batch/BELLPACK_a2 $path/${Matrices[$i]} $row_blk_size

          #           echo "sbatch batch/CUTLASS_BELLPACK_a5 $path/${Matrices[$i]} $row_blk_size ${Taus[$i]}"
#           sbatch batch/CUTLASS_BELLPACK_a5 $path/${Matrices[$i]} $row_blk_size ${Taus[$i]}

#           echo "sbatch batch/CUTLASS_BELLPACK_a2 $path/${Matrices[$i]} $row_blk_size"
          sbatch batch/CUTLASS_BELLPACK_a2 $path/${Matrices[$i]} $row_blk_size
          fi

          if [ $row_blk_size -eq 32 ] && [ $row_blk_size -eq 32 ]
          then
#               echo "sbatch batch/CSR $path/${Matrices[$i]}"
            sbatch batch/CSR $path/${Matrices[$i]}

#               echo "sbatch batch/CSR $path/${Matrices[$i]}"
            sbatch batch/DENSE $path/${Matrices[$i]}
          fi
        else
          echo "ERROR: $path/${Matrices[$i]} does not exists."
        fi

    done

done
