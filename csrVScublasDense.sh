#!/bin/bash

path=../SPARTA_datasets/uniform_sets

Matrices=('uniform_1024x1024_0001.el'  'uniform_1024x1024_0050.el'  'uniform_1024x1024_2000.el' 'uniform_8192x8192_0001.el'  'uniform_8192x8192_0050.el' 'uniform_1024x1024_0002.el'  'uniform_1024x1024_0100.el'  'uniform_32768x32768_0001.el'  'uniform_8192x8192_0002.el'  'uniform_8192x8192_0100.el' 'uniform_1024x1024_0005.el'  'uniform_1024x1024_0200.el'  'uniform_32768x32768_0002.el'  'uniform_8192x8192_0005.el'  'uniform_8192x8192_0200.el' 'uniform_1024x1024_0010.el'  'uniform_1024x1024_0500.el'  'uniform_32768x32768_0005.el'  'uniform_8192x8192_0010.el' 'uniform_1024x1024_0020.el'  'uniform_1024x1024_1000.el'  'uniform_32768x32768_0010.el'  'uniform_8192x8192_0020.el')


for i in ${!Matrices[@]}
do
#  echo "Matrix ${Matrices[$i]}"

  for c in 2048 8192
  do

    if [ -f "$path/${Matrices[$i]}" ]
    then

      #  echo "sbatch batch/CSR $path/${Matrices[$i]} $c"
      sbatch batch/CSR $path/${Matrices[$i]} 0 $c

      #  echo "sbatch batch/DENSE $path/${Matrices[$i]} $c"
      sbatch batch/DENSE $path/${Matrices[$i]} 0 $c

    else
      echo "ERROR: $path/${Matrices[$i]} does not exists."
    fi

  done
done
