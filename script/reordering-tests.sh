#!/bin/bash

OPTS="-r 3 -v -1 -i 4"

: > "test_reordering_synth.txt"

m_shapes=(1024 2048);
n_shapes=(1024 2048);

P_value=(23 41 64 128 200 256);
p_value=(8 16 32 64 128 200);
b_value=(0.01 0.1 0.2 0.3 0.5 0.7 0.9);
q_value=(0.01 0.02 0.05 0.1 0.3 0.5 0.8 1.0);
e_value=(0.5 0.6 0.65 0.7 0.75 0.8);
s_value=(1 3);

for m in ${m_shapes[@]}; do
  for n in ${n_shapes[@]}; do
    for e in ${e_value[@]}; do
      for p in ${p_value[@]}; do
	for P in ${P_value[@]}; do
          for b in ${b_value[@]}; do
            for q in ${q_value[@]}; do
	      for s in ${s_value[@]}; do
		  echo $m $n $e $p $P $b $q
	          ./programs/general/test_AHA $OPTS -m $m -n $n -e $e -p $p -P $P -b $b -q $q >> "test_reordering_synth.txt"
              done
            done
	  done
        done
      done
    done
  done
done



