#!/bin/bash


MATRIX_PATH=$1
mkdir $MATRIX_PATH


nmin=20000
nmax=100000
increment=5000

dmin=10000
dmax=50

array=()

for ((i=nmin; i<=nmax; i+=increment)); do
	nmin_tmp=$i
	nmax_tmp=$((i+increment))
	nzmin=$((nmin_tmp*nmin_tmp/dmin))
	nzmax=$((nmax_tmp*nmax_tmp/dmax))
	ssget_line="[ @cols -ge ${nmin_tmp} ] && [ @cols -le ${nmax_tmp} ] && [ @rows -ge ${nmin_tmp} ] && [ @rows -le ${nmax_tmp} ] && [ @nonzeros -ge ${nzmin} ] && [ @nonzeros -le ${nzmax} ]"
	tmp_array=($(ssget -s "${ssget_line}"))
	array+=("${tmp_array[@]}")
done

for id in "${array[@]}";do
	echo $(ssget -i $id -p name)
	filename=$(ssget -f -i $id -e -t MM)
	cp $filename $MATRIX_PATH
done
