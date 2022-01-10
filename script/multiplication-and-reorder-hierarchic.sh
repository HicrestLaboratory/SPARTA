#!/bin/bash

OPTS="-r 3 -v -1 -i 4 -R saad_blocks -s 1 -M 0 -F jaccard"


DATE=$(date +"%m-%d-%Y");
RESULTS="results/test_cublas_reordering-non-hierarchic-${DATE}.txt"

:>${RESULTS};

m_shapes=(8192);
k_shapes=(8192);
n_shapes=(128);

e_value=(0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9);

P_value=(64);

b_value=(0.01 0.05 0.1 0.2 0.3 0.4 0.5);
q_value=(0.005 0.01 0.02 0.05 0.1 0.2 0.5);

l_value=(0);

for m in ${m_shapes[@]}; do
  for k in ${k_shapes[@]}; do
    for e in ${e_value[@]}; do
	for P in ${P_value[@]}; do
          for b in ${b_value[@]}; do
            for q in ${q_value[@]}; do
	      for n in ${n_shapes[@]}; do
		  for l in ${l_value[@]}; do 
		    echo "m,n,k " $m $n $k "e,p,b,q,l " $e $P $b $q $l
	            ./programs/cuda/test_cublas_reorder ${OPTS} -m $m -k $k -n $n -e $e -p $P -P $P -b $b -q $q -l $l >> ${RESULTS}
	 	  done
	      done
            done
	  done
        done
    done
  done
done



