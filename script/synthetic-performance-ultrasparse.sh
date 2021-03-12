#!/bin/bash

OPTS="-w 1 -r 3 -v -1 -i 4"


RESULTS = "results/test_complete_cublas_results_ultrasparse.txt
: > $(RESULTS)



Mshapes=(1024 2048 4096 8192 16384);
Nshapes=(1024 2048 4096 8192 16384);
Kshapes=(1024 2048 4096 8192 16384);

Pvalue=(32 64 128 256 512 1024);
Bvalue=(0.05 0.1 0.15 0.2 0.25 0.3 0.4);
Qvalue=(0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.05 0.1);

for m in ${Mshapes[@]}; do
  for k in ${Kshapes[@]}; do
    for n in ${Nshapes[@]}; do
      for p in ${Pvalue[@]}; do
        for b in ${Bvalue[@]}; do
          for q in ${Qvalue[@]}; do
            	./programs/cuda/test_cublas_VBS $OPTS -m $m -n $n -k $k -p $p -b $b -q $q >> $(RESULTS)
          done
        done
      done
    done
  done
done



