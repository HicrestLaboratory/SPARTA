#!/bin/bash

path=../SPARTA_datasets/suitsparse_collection_3_no3-5

Matrices=('ckt11752_dc_1.mtx' 'email-Enron.mtx' 'FEM_3D_thermal1.mtx' 'g7jac100sc.mtx' 'k3plates.mtx' 'lhr17.mtx' 'TEM27623.mtx')
Taus=(0.2  0.05 0.01 0.01 0.9  0.8  0.5 )

row_blk_size=128
col_blk_size=512

for i in ${!Matrices[@]}; do
  echo "With $row_blk_size x $col_blk_size, matrix ${Matrices[$i]} has tau ${Taus[$i]}"
  echo "sbatch batch/VBR $path/${Matrices[$i]} 128 512 ${Taus[$i]}"
  sbatch batch/VBR_rect $path/${Matrices[$i]} $row_blk_size $col_blk_size ${Taus[$i]}

  echo "sbatch batch/CSR $path/${Matrices[$i]}"
  sbatch batch/CSR $path/${Matrices[$i]}
done

echo "-----------------------------------------------------------------------"

Matrices=('ckt11752_dc_1.mtx' 'email-Enron.mtx' 'FEM_3D_thermal1.mtx' 'g7jac100sc.mtx' 'k3plates.mtx' 'lhr17.mtx' 'poisson3Da.mtx' 'TEM27623.mtx' 'vsp_south31_slptsk.mtx')
Taus=(0.9  0.05 0.1  0.01 0.4  0.9  0.01 0.9  0.01)

row_blk_size=512
col_blk_size=512

for i in ${!Matrices[@]}; do
  echo "With $row_blk_size x $col_blk_size, matrix ${Matrices[$i]} has tau ${Taus[$i]}"
  echo "sbatch batch/VBR $path/${Matrices[$i]} 128 512 ${Taus[$i]}"
  sbatch batch/VBR_rect $path/${Matrices[$i]} $row_blk_size $col_blk_size ${Taus[$i]}

  echo "sbatch batch/CSR $path/${Matrices[$i]}"
  sbatch batch/CSR $path/${Matrices[$i]}
done

echo "-----------------------------------------------------------------------"

path=../SPARTA_datasets/suitsparse_collection_4

Matrices=('ca-CondMat.mtx' 'cz20468.mtx' 'epb2.mtx' 'jan99jac100sc.mtx' 'lowThrust_5.mtx' 'rajat22.mtx' 'vsp_barth5_1Ksep_50in_5Kout.mtx')
Taus=(0.8  0.4  0.6  0.05 0.8  0.2  0.9)

row_blk_size=128
col_blk_size=512

for i in ${!Matrices[@]}; do
  echo "With $row_blk_size x $col_blk_size, matrix ${Matrices[$i]} has tau ${Taus[$i]}"
  echo "sbatch batch/VBR $path/${Matrices[$i]} 128 512 ${Taus[$i]}"
  sbatch batch/VBR_rect $path/${Matrices[$i]} $row_blk_size $col_blk_size ${Taus[$i]}

  echo "sbatch batch/CSR $path/${Matrices[$i]}"
  sbatch batch/CSR $path/${Matrices[$i]}
done

echo "-----------------------------------------------------------------------"

Matrices=('ca-CondMat.mtx' 'cz20468.mtx' 'epb2.mtx' 'jan99jac100sc.mtx' 'lowThrust_5.mtx' 'rajat22.mtx' 'vsp_barth5_1Ksep_50in_5Kout.mtx')
Taus=(0.8  0.01 0.5  0.1  0.9  0.3  0.9)

row_blk_size=512
col_blk_size=512

for i in ${!Matrices[@]}; do
  echo "With $row_blk_size x $col_blk_size, matrix ${Matrices[$i]} has tau ${Taus[$i]}"
  echo "sbatch batch/VBR $path/${Matrices[$i]} 128 512 ${Taus[$i]}"
  sbatch batch/VBR_rect $path/${Matrices[$i]} $row_blk_size $col_blk_size ${Taus[$i]}

  echo "sbatch batch/CSR $path/${Matrices[$i]}"
  sbatch batch/CSR $path/${Matrices[$i]}
done
