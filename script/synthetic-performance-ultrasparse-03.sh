#!/bin/bash

OPTS="-w 1 -r 10 -v -1 -i 4"



Mshapes=(2048 4096 8192);
Nshapes=(2048 4096 8192 16384);
Kshapes=(2048 4096 8192);

Pvalue=(32 64 128);
Bvalue=(0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5);
Qvalue=(0.03);
RESULTS="results/test_cublas_results_ultrasparse_19_10_q"-$Qvalue
:>"${RESULTS}";

for m in ${Mshapes[@]}; do
  for k in ${Kshapes[@]}; do
    for n in ${Nshapes[@]}; do
      for p in ${Pvalue[@]}; do
        for b in ${Bvalue[@]}; do
          for q in ${Qvalue[@]}; do
          	./programs/cuda/test_cublas_VBS ${OPTS} -m $m -n $n -k $k -p $p -b $b -q\ $q >> "${RESULTS}"
          done
        done
      done
    done
  done
done



