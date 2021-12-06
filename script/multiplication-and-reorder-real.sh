#!/bin/bash

OPTS="-r 3 -v -1 -i 2 -R saad_blocks -s 1 -M 1 -a 1 -F jaccard"
DATA={1}

DATE=$(date +"%m-%d-%Y");
RESULTS="results/test_cublas_reordering-${DATE}.txt"

:>${RESULTS};

n_shapes=(2048 4096 8192)

e_value=(0.1 0.2 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95);

P_value=(32 64 128);
l_value=(0 1 2);

for f in ${1}*; do
  for e in ${e_value[@]}; do
    for P in ${P_value[@]}; do
      for n in ${n_shapes[@]}; do
        for l in ${l_value[@]}; do 
	      echo $f $P
	      ./programs/cuda/test_cublas_reorder ${OPTS} -f $f -n $n -e $e -P $P -l $l >> ${RESULTS}
	done
      done
    done
  done
done



