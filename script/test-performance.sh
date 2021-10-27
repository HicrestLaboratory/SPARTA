#!/bin/bash

OPTS="-w 1 -r 5 -v -1 -i 4"

RESULTS="results/small_test.txt";

>"${RESULTS}";

Mshapes=(8192);
Nshapes=(8192);
Kshapes=(8192);

Pvalue=(64);
Bvalue=(0.1);
Qvalue=(0.1);

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



