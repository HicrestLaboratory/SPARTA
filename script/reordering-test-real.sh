#!/bin/bash

OPTS="-i 2 -v -1"
​
: > "test_inputs.txt"

P_value=(5 8 12 16 32 64 128 256);
e_value=(0.5 0.6 0.8 0.9 0.99);
s_value=(1 3);
​
for f in "data/"*; do
	for e in ${e_value[@]}; do
		for P in ${P_value[@]}; do
	      		for s in ${s_value[@]}; do
				echo $f $e $P $s
	          		./programs/general/test_AHA $OPTS -e $e -P $P -s $s -f $f >> "test_inputs.txt"
			done
		done
	done
done
