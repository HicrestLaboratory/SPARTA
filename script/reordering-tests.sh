#!/bin/bash

OPTS="-r 5 -v -1 -i 4"

: > "test_reordering.txt"

m_shapes=(1024 2048 4096 8192);
n_shapes=(1024 2048 4096 8192);

P_value=(32 35 50 64 100 256 512);
p_value=(64 128 256 512 1024);
b_value=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0);
q_value=(0.001 0.01 0.05 0.1 0.15 0.20 0.25 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0);
e_value= (0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99);

for m in ${m_shapes[@]}; do
  for n in ${n_shapes[@]}; do
    for e in ${e_value[@]}; do
      for p in ${p_value[@]}; do
	for p in ${P_value[@]}; do
          for b in ${b_value[@]}; do
            for q in ${q_value[@]}; do
              ./programs/general/test_AHA $OPTS -m $m -n $n -e $e -p $p -P $P -b $b -q $q >> "test_reordering.txt"
            done
	  done
        done
      done
    done
  done
done



