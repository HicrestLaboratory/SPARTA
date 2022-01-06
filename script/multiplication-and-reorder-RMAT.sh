#!/bin/bash

OPTS="-r 4 -v -1 -i 2 -R saad_blocks -s 1 -a -1 -F jaccard -M 1"
DATA={1}

DATE=$(date +"%m-%d-%Y");
RESULTS="results/test_cublas_reordering-RMAT-${DATE}.txt"

:>${RESULTS};

e_value=(0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9);

P_value=(64 128 256);
l_value=(-1 0 2 4);

for f in ${1}*; do
  for e in ${e_value[@]}; do
    for P in ${P_value[@]}; do
      for l in ${l_value[@]}; do 
		echo $f $e $P $l $(date)
	      	./programs/cuda/test_cublas_reorder_optimized ${OPTS} -f $f -e $e -P $P -l $l >> ${RESULTS}
      done
    done
  done
done



