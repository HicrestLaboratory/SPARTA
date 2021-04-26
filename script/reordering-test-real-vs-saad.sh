#!/bin/bash

OPTS="-i 2 -v -1"

RESULTS="results/reordering_vs_saad.txt";

>"${RESULTS}";

P_value=(16 32 64 128 200 256 300 512);
e_value=(0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9);
s_value=(1 3);

for f in data/real-graphs/*.txt; do
	for e in ${e_value[@]}; do
		for P in ${P_value[@]}; do
	      		for s in ${s_value[@]}; do
				echo $f $e $P $s
	          		./programs/general/test_saad $OPTS -e $e -P $P -s $s -f $f -R saad_blocks >> "${RESULTS}"
				./programs/general/test_saad $OPTS -e $e -P $P -s $s -f $f -R saad >> "${RESULTS}"
			done
		done
	done
done
