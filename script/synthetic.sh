#!/bin/bash

OPTS="-w 1 -r 10 -v -1 -i 4"

Output_file="test_cublas_results.txt"

Mshapes=(1024 2048 4096 8192 16384);
Nshapes=(1024 2048 4096 8192 16384);
Kshapes=(1024 2048 4096 8192 16384);

Pvalue=(128 256 512 1024);
Bvalue=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0);
Qvalue=(0.001 0.01 0.05 0.1 0.15 0.20 0.25 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0);

for m in ${Mshapes[@]}; do
  for k in ${Kshapes[@]}; do
    for n in ${Nshapes[@]}; do
#      #echo "$m $k $n"; sleep 0.2 
      for p in ${Pvalue[@]}; do
        for b in ${Bvalue[@]}; do
          for q in ${Qvalue[@]}; do
            ./programs/cuda/test_cublas_VBS $OPTS -m $m -n $n -k $k -p $p -b $b -q $q; sleep 0.5; >> ${Output_file}
          done
        done
      done
    done
  done
done




