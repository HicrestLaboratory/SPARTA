#!/bin/bash

OPTS="-r 5 -v -1 -i 4 -R saad_blocks"

RESULTS="results/test_reordering_blocked_synth_12_11.txt"
:>${RESULTS};

m_shapes=(2048);
k_shapes=(2048);

e_value=(0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.999);

P_value=(32 64 128);

b_value=(0.1 0.2 0.3 0.4 0.5);
q_value=(0.01 0.02 0.05 0.1 0.2 0.5);

F_value=(scalar jaccard);
M_value=(0 1);
l_value=(-1 0.25 0.5 1 2);


for m in ${m_shapes[@]}; do
  for k in ${k_shapes[@]}; do
    for e in ${e_value[@]}; do
	for P in ${P_value[@]}; do
          for b in ${b_value[@]}; do
            for q in ${q_value[@]}; do
	      for M in ${M_value[@]}; do
		for F in ${F_value[@]}; do
		  for l in ${l_value[@]}; do 
		    echo $m $n $e $p $P $b $q $F $M $l
	            ./programs/general/test_saad ${OPTS} -m $m -k $k -e $e -p $P -P $P -b $b -q $q -F $F -M $M -l $l>> ${RESULTS}
	 	  done
                done
	      done
            done
	  done
        done
    done
  done
done



