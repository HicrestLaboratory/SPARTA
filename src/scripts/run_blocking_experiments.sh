export MATRICES_PATH=$1
export RESULTS_PATH=$2
export PROGRAM=$3


TAUs=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
USE_PATTERN=(0 1)
USE_STRUCTURED=(0 1)
USE_GROUP=(0 1)
BLOCK_SIZEs=(32 64 128)

for fullpath in ${MATRICES_PATH}/*.el; do
	MATRIX_FILE=$(basename -- "${fullpath}")
	MATRIX_NAME="${MATRIX_FILE%.*}"
	echo "============= processing matrix ${MATRIX_NAME}"
	MATRIX_FOLDER=${RESULTS_PATH}/${MATRIX_NAME}
	mkdir ${MATRIX_FOLDER}
	for b in ${BLOCK_SIZEs[@]}; do
		for S in ${USE_STRUCTURED[@]}; do
			for t in ${TAUs[@]}; do
				for p in ${USE_PATTERN[@]}; do
					for g in ${USE_GROUP[@]}; do
						EXP_NAME="blocking_G_${MATRIX_NAME}_b_${b}_t_${t}_S_${S}_p_${p}_g_${g}"
						OUTFILE=${MATRIX_FOLDER}/${EXP_NAME}.txt
						./${PROGRAM} -f ${fullpath} -b ${b} -t ${t} -S ${S} -p ${p} -g ${g} -v 1 -o ${OUTFILE} -P 1 -n ${EXP_NAME}
					done
				done
			done
		done
	done


done
