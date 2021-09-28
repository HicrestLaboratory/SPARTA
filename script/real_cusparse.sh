#!/bin/bash

OPTS="-w 1 -r 10 -v -1 -i 2 -a 5"

RESULTS = "results/real_cusparse_results.txt"
: > ${RESULTS}

n_shapes=(1024 2048 4096 8192 16384);

for f in "data/*"; do
	for n in ${n_shapes[@]}; do
            	./programs/cuda/test_cublas_VBS $OPTS -f $f -n $n>> ${RESULTS}
	done
done



