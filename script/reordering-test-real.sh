#!/bin/bash

OPTS="-i 2 -v -1"

P_value=(64 128 200 256 300 512);
e_value=(0.5 0.6 0.65 0.7 0.75 0.8 );
s_value=(1 3);
​
RESULTS = "results/real_reordering_results.txt"
: > "$(RESULTS)"

for f in "data/*"; do
	for e in ${e_value[@]}; do
		for P in ${P_value[@]}; do
	      		for s in ${s_value[@]}; do
				echo $f $e $P $s
	          		./programs/general/test_AHA $OPTS -e $e -P $P -s $s -f $f >> $(RESULTS)
			done
		done
	done
done
