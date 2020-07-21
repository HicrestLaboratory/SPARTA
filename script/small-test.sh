#!/bin/bash

OPTS="-w 1 -r 3 -v -1 -i 4"

: > test_small_cublas_results.txt

Mshapes=(1024 2048);
Nshapes=(1024 2048);
Kshapes=(1024 2048);

Pvalue=(128 512)
Bvalue=(0.1 0.5 1);
Qvalue=(0.1 0.5 1.0);

for m in ${Mshapes[@]}; do
  for k in ${Kshapes[@]}; do
    for n in ${Nshapes[@]}; do
#      #echo "$m $k $n"; sleep 0.2 
      for p in ${Pvalue[@]}; do
        for b in ${Bvalue[@]}; do
          for q in ${Qvalue[@]}; do
            ./programs/cuda/test_cublas_VBS $OPTS -m $m -n $n -k $k -p $p -b $b -q $q; sleep 0.5; 1 >> test_small_cublas_results.txt
          done
        done
      done
    done
  done
done
