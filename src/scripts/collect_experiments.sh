SOURCES=("$@")
outputfile=${SOURCES[0]}
unset SOURCES[0]


first=0


first_folder=${SOURCES[1]}
for csv_file in ${first_folder}/*/*.*; do
	sed -n '1 p' "${csv_file}" > ${outputfile}
	break
done

for sou_folder in ${SOURCES[@]}; do
	for csv_file in ${sou_folder}/*/*.*; do
		sed -n '2 p' "${csv_file}" >> "${outputfile}"
	done
done




