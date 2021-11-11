#!/bin/bash
OPTS="-i 1 -v -1 -w 1 -r 5";
RESULTS="results/cublas_cusparse_comparison.txt";

>"${RESULTS}";

m_value=(128 256 512 1024 2048 4096 8192);
n_value=(128 256 512 1024 2048 4096 8192);
k_value=(128 256 512 1024 2048 4096 8192);
q_value=(0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.15 0.2 0.25)  

for m in ${m_value[@]}; do
	for n in ${n_value[@]}; do
		for k in ${k_value[@]}; do
	      		for q in ${q_value[@]}; do
				echo $m $n $k $q
	          		/home/clusterusers/pasyloslabini/SPARTA/programs/cuda/test_cublas_cusparse_comparison ${OPTS} -m $m -n $n -k $k -q $q>> "${RESULTS}"
			done
		done
	done
done
