#!/bin/bash

OPTS="-w 1 -r 10 -v -1 -i 4"


Mshapes=(1024 4096);
Nshapes=(1024 4096);
Kshapes=(1024 4096);

Pvalue=(256 512);
Bvalue=(0.1 0.5 0.9);
Qvalue=(0.05 0.1 0.4 0.8);

for m in ${Mshapes[@]}; do
  for k in ${Kshapes[@]}; do
    for n in ${Nshapes[@]}; do
#      #echo "$m $k $n"; sleep 0.2 
      for p in ${Pvalue[@]}; do
        for b in ${Bvalue[@]}; do
          for q in ${Qvalue[@]}; do
            ./cuda_test $OPTS -m $m -n $n -k $k -p $p -b $b -q $q sleep 0.5;
          done
        done
      done
    done
  done
done




