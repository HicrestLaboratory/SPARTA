#!/bin/bash

OPTS="-i 2 -v -1 -R saad-blocks -s 1"

P_value=(64 128 200 256 300 512);
e_value=(0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9);
F_value=(scalar jaccard);
<200b>
RESULTS="results/real_reordering_results_25_10.txt";

for f in "data/*"; do
        for e in ${e_value[@]}; do
                for P in ${P_value[@]}; do
                        for F in ${F_value[@]}; do
                                echo $f $e $P $F
                                ./programs/general/test_saad ${OPTS} -e $e -P $P -F $F>> ${RESULTS}
                        done
                done
        done
done
~                                                                                                                                                                                            ~                                                                                                                                                                                            ~                                                                                                                                                                                            ~                                                                                                                                                                                            ~                                                                                                                                                                                            ~                                                                                                                                                                                            ~                                            