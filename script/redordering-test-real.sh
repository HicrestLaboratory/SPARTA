#!/bin/bash
​  
arr=(../data/*)
for f in arr; do
	echo f

OPTS="-i 2 -f" + $(path) + $(filename) + "-v -1"
​
: > "test_inputs.txt"
​
P_value=(8 16 32 64 128 256);
e_value=(0.5 0.6 0.8 0.9 0.99);
s_value=(1 3);
​
for e in ${e_value[@]}; do
	for P in ${P_value[@]}; do
	      	for s in ${s_value[@]}; do
			echo $e $P $s
	          	./../programs/general/test_AHA $OPTS -e $e -P $P -s $s >> "test_inputs.txt"
		done
	done
done