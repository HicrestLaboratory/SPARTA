#!/bin/bash

OPTS="-w 1 -r 10 -v -1 -i 2 -a 5"

RESULTS = "results/real_cusparse_results.txt"
: > $(RESULTS)

for f in "data/*"; do
            	./programs/cuda/test_cublas_VBS $OPTS -f $f >> $(RESULTS)
done



