#!/bin/bash

OPTS="-i 2 -v -1"

P_value=(64 128 200 256 300 512);
e_value=(0.7 0.8 0.9);
s_value=(1 3);
​
RESULTS="results/compression_results.txt";
: > $(RESULTS);

for f in "data/*"; do
	for e in ${e_value[@]}; do
		for P in ${P_value[@]}; do
	      		for s in ${s_value[@]}; do
				echo $f $e $P $s
	          		./programs/general/test_M1vsM2 $OPTS -e $e -P $P -s $s -f $f >> $(RESULTS)
			done
		done
	done
done
