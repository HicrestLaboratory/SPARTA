#!/bin/bash

path=../SPARTA_datasets/suitsparse_collection_3

Matrices=('vsp_south31_slptsk.mtx' 'ckt11752_dc_1.mtx' 'email-Enron.mtx' 'TEM27623.mtx' 'g7jac100sc.mtx' 'lhr17.mtx' 'FEM_3D_thermal1.mtx' 'poisson3Da.mtx' 'k3plates.mtx')
Taus=( 0.2  0.2  0.2  0.2  0.2 0.2  0.2  0.2  0.2 )

for row_blk_size in 32 64 128 256 512 1024
do
  for col_blk_size in 32 64 128 256 512 1024
  do

    echo "row x col = $row_blk_size x $col_blk_size"
    if [[ "$row_blk_size" -eq 32 ]] && [[ "$col_blk_size" -eq 64 ]]
    then
      Taus=( 0.4  0.3  0.2  0.9  0.05 0.3  0.4  0.2  0.6 )
    fi
    if [[ "$row_blk_size" -eq 32 ]] && [[ "$col_blk_size" -eq 128 ]]
    then
     Taus=( 0.4  0.3  0.2  0.9  0.05 0.4  0.3  0.2  0.4 )
    fi
    if [[ "$row_blk_size" -eq 32 ]] && [[ "$col_blk_size" -eq 256 ]]
    then
     Taus=( 0.3  0.3  0.3  0.2  0.05 0.1  0.05 0.2  0.4 )
    fi
    if [[ "$row_blk_size" -eq 32 ]] && [[ "$col_blk_size" -eq 512 ]]
    then
     Taus=( 0.2  0.3  0.05 0.01 0.01 0.05 0.05 0.1  0.4 )
    fi
    if [[ "$row_blk_size" -eq 32 ]] && [[ "$col_blk_size" -eq 1024 ]]
    then
     Taus=( 0.2  0.2  0.2  0.01 0.01 0.01 0.01 0.01 0.2 )
    fi
    if [[ "$row_blk_size" -eq 64 ]] && [[ "$col_blk_size" -eq 64 ]]
    then
     Taus=( 0.2  0.9  0.3  0.9  0.05 0.4  0.9  0.2  0.9 )
    fi
    if [[ "$row_blk_size" -eq 64 ]] && [[ "$col_blk_size" -eq 128 ]]
    then
     Taus=( 0.4  0.9  0.2  0.9  0.05 0.3  0.3  0.2  0.8 )
    fi
    if [[ "$row_blk_size" -eq 64 ]] && [[ "$col_blk_size" -eq 256 ]]
    then
     Taus=( 0.3  0.3  0.1  0.9  0.01 0.1  0.05 0.2  0.6 )
    fi
    if [[ "$row_blk_size" -eq 64 ]] && [[ "$col_blk_size" -eq 512 ]]
    then
     Taus=( 0.2  0.3  0.1  0.01 0.01 0.1  0.05 0.1  0.9 )
    fi
    if [[ "$row_blk_size" -eq 64 ]] && [[ "$col_blk_size" -eq 1024 ]]
    then
     Taus=( 0.3  0.2  0.01 0.01 0.01 0.05 0.01 0.01 0.4 )
    fi
    if [[ "$row_blk_size" -eq 128 ]] && [[ "$col_blk_size" -eq 64 ]]
    then
     Taus=( 0.3  0.9  0.05 0.9  0.05 0.9  0.9  0.1  0.9 )
    fi
    if [[ "$row_blk_size" -eq 128 ]] && [[ "$col_blk_size" -eq 128 ]]
    then
     Taus=( 0.2  0.9  0.4  0.9  0.01 0.9  0.3  0.1  0.8 )
    fi
    if [[ "$row_blk_size" -eq 128 ]] && [[ "$col_blk_size" -eq 256 ]]
    then
     Taus=( 0.1  0.9  0.05 0.9  0.05 0.9  0.3  0.2  0.8 )
    fi
    if [[ "$row_blk_size" -eq 128 ]] && [[ "$col_blk_size" -eq 512 ]]
    then
     Taus=( 0.1  0.2  0.05 0.9  0.01 0.8  0.01 0.1  0.9 )
    fi
    if [[ "$row_blk_size" -eq 128 ]] && [[ "$col_blk_size" -eq 1024 ]]
    then
     Taus=( 0.1  0.2  0.1  0.01 0.01 0.01 0.01 0.01 0.4 )
    fi
    if [[ "$row_blk_size" -eq 256 ]] && [[ "$col_blk_size" -eq 64 ]]
    then
     Taus=( 0.4 0.9 0.1 0.9 0.1 0.9 0.8 0.1 0.9 )
    fi
    if [[ "$row_blk_size" -eq 256 ]] && [[ "$col_blk_size" -eq 128 ]]
    then
     Taus=( 0.05 0.9  0.05 0.9  0.01 0.9  0.9  0.1  0.8 )
    fi
    if [[ "$row_blk_size" -eq 256 ]] && [[ "$col_blk_size" -eq 256 ]]
    then
     Taus=( 0.05 0.9  0.05 0.9  0.01 0.8  0.2  0.05 0.4 )
    fi
    if [[ "$row_blk_size" -eq 256 ]] && [[ "$col_blk_size" -eq 512 ]]
    then
     Taus=( 0.3  0.3  0.05 0.9  0.01 0.8  0.1  0.05 0.4 )
    fi
    if [[ "$row_blk_size" -eq 256 ]] && [[ "$col_blk_size" -eq 1024 ]]
    then
     Taus=( 0.01 0.2  0.05 0.3  0.01 0.01 0.01 0.01 0.4 )
    fi
    if [[ "$row_blk_size" -eq 512 ]] && [[ "$col_blk_size" -eq 64 ]]
    then
     Taus=( 0.01 0.9  0.2  0.9  0.1  0.9  0.8  0.05 0.9 )
    fi
    if [[ "$row_blk_size" -eq 512 ]] && [[ "$col_blk_size" -eq 128 ]]
    then
     Taus=( 0.05 0.9  0.2  0.9  0.01 0.8  0.9  0.05 0.8 )
    fi
    if [[ "$row_blk_size" -eq 512 ]] && [[ "$col_blk_size" -eq 256 ]]
    then
     Taus=( 0.1  0.9  0.05 0.9  0.01 0.9  0.3  0.1  0.8 )
    fi
    if [[ "$row_blk_size" -eq 512 ]] && [[ "$col_blk_size" -eq 512 ]]
    then
     Taus=( 0.01 0.9  0.05 0.9  0.01 0.9  0.1  0.01 0.4 )
    fi
    if [[ "$row_blk_size" -eq 512 ]] && [[ "$col_blk_size" -eq 1024 ]]
    then
     Taus=( 0.01 0.05 0.01 0.8  0.2  0.2  0.05 0.01 0.05 )
    fi
    if [[ "$row_blk_size" -eq 1024 ]] && [[ "$col_blk_size" -eq 64 ]]
    then
     Taus=( 0.3 0.8 0.3 0.9 0.9 0.9 0.8 0.1 0.9 )
    fi
    if [[ "$row_blk_size" -eq 1024 ]] && [[ "$col_blk_size" -eq 128 ]]
    then
     Taus=( 0.01 0.9  0.2  0.9  0.4  0.8  0.9  0.05 0.8 )
    fi
    if [[ "$row_blk_size" -eq 1024 ]] && [[ "$col_blk_size" -eq 256 ]]
    then
     Taus=( 0.01 0.9  0.1  0.9  0.01 0.8  0.01 0.05 0.8 )
    fi
    if [[ "$row_blk_size" -eq 1024 ]] && [[ "$col_blk_size" -eq 512 ]]
    then
     Taus=( 0.05 0.9  0.3  0.9  0.01 0.8  0.3  0.01 0.8 )
    fi
    if [[ "$row_blk_size" -eq 1024 ]] && [[ "$col_blk_size" -eq 1024 ]]
    then
     Taus=( 0.9  0.8  0.05 0.8  0.2  0.8  0.05 0.01 0.01)
    fi

    for i in ${!Matrices[@]}
    do
        echo "With $row_blk_size x $col_blk_size, matrix ${Matrices[$i]} has tau ${Taus[$i]}"

        if [ -f "$path/${Matrices[$i]}" ]
        then

#           echo "sbatch batch/VBR_rect_a5 $path/${Matrices[$i]} $row_blk_size $col_blk_size ${Taus[$i]}"
          sbatch batch/VBR_rect_a5 $path/${Matrices[$i]} $row_blk_size $col_blk_size ${Taus[$i]}

#           echo "sbatch batch/VBR_rect_a2 $path/${Matrices[$i]} $row_blk_size $col_blk_size"
          sbatch batch/VBR_rect_a2 $path/${Matrices[$i]} $row_blk_size $col_blk_size

          #           echo "sbatch batch/VBR_rect_a5 $path/${Matrices[$i]} $row_blk_size $col_blk_size ${Taus[$i]}"
          sbatch batch/VBR_batch_a5 $path/${Matrices[$i]} $row_blk_size $col_blk_size ${Taus[$i]}

#           echo "sbatch batch/VBR_rect_a2 $path/${Matrices[$i]} $row_blk_size $col_blk_size"
          sbatch batch/VBR_batch_a2 $path/${Matrices[$i]} $row_blk_size $col_blk_size

          if [ $row_blk_size == $col_blk_size ]
          then
#           echo "sbatch batch/BELLPACK_a5 $path/${Matrices[$i]} $row_blk_size ${Taus[$i]}"
          sbatch batch/BELLPACK_a5 $path/${Matrices[$i]} $row_blk_size ${Taus[$i]}

          #           echo "sbatch batch/BELLPACK_a2 $path/${Matrices[$i]} $row_blk_size"
          sbatch batch/BELLPACK_a2 $path/${Matrices[$i]} $row_blk_size

          #           echo "sbatch batch/CUTLASS_BELLPACK_a5 $path/${Matrices[$i]} $row_blk_size ${Taus[$i]}"
          sbatch batch/CUTLASS_BELLPACK_a5 $path/${Matrices[$i]} $row_blk_size ${Taus[$i]}

#           echo "sbatch batch/CUTLASS_BELLPACK_a2 $path/${Matrices[$i]} $row_blk_size"
          sbatch batch/CUTLASS_BELLPACK_a2 $path/${Matrices[$i]} $row_blk_size
          fi

          if [ $row_blk_size -eq 32 ] && [ $col_blk_size -eq 32 ]
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
done
