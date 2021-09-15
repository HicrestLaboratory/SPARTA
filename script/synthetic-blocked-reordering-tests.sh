#!/bin/bash

OPTS="-r 3 -v -1 -i 4 -R saad-blocks -N ReorderingBlockedSynth"

RESULTS = "results/test_reordering_blocked_synth.txt"
: > $(RESULTS)

m_shapes=(1024 2048);
k_shapes=(1024 2048);

e_value=(0.5 0.6 0.65 0.7 0.75 0.8);

P_value=(16 23 41 64 100 128);
p_value=(16 32 64 128);

b_value=(0.01 0.1 0.2 0.3 0.5 0.7 0.9);
q_value=(0.005 0.01 0.02 0.05 0.1 0.3);

F_value=(jaccard scalar);
s_value=(1 3);

for m in ${m_shapes[@]}; do
  for k in ${k_shapes[@]}; do
    for e in ${e_value[@]}; do
      for p in ${p_value[@]}; do
	for P in ${P_value[@]}; do
          for b in ${b_value[@]}; do
            for q in ${q_value[@]}; do
	      for s in ${s_value[@]}; do
		for F in ${F_value[@]}; do
		  echo $m $n $e $p $P $b $q
	          ./programs/general/test_AHA $OPTS -m $m -k $k -e $e -p $p -P $P -b $b -q $q -F $F -s $s>> $(RESULTS)
                done
	      done
            done
	  done
        done
      done
    done
  done
done



