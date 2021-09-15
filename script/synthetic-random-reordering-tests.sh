#!/bin/bash

OPTS="-r 3 -v -1 -i 1 -s 1 -R saad-blocks -N ReorderingSynthRand"

RESULTS = "results/test_reordering_blocked_synth.txt"
: > ${RESULTS};

m_shapes=(1024 2048);
k_shapes=(1024 2048);

P_value=(16 23 41 64 100 128);
q_value=(0.001 0.005 0.01 0.02 0.05 0.1);
e_value=(0.5 0.6 0.65 0.7 0.75 0.8);
F_value=(jaccard scalar);

for m in ${m_shapes[@]}; do
  for k in ${k_shapes[@]}; do
    for e in ${e_value[@]}; do
      for P in ${P_value[@]}; do
        for q in ${q_value[@]}; do
	  for F in ${F_value[@]}; do
	    echo $m $n $e $p $P $b $q
	    ./programs/general/test_AHA $OPTS -m $m -k $k -e $e -P $P -q $q -F $F>> ${RESULTS}
          done
        done
      done
    done
  done
done



