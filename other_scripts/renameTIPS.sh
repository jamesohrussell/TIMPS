IFS='_'
for fil in TIMPS*.nc
do
  echo $fil
  read -a strarr <<< "$fil"
  mv "$fil" "${strarr[0]}_`printf "%07d" ${strarr[1]}`_${strarr[2]}_${strarr[3]}_${strarr[4]}"
done
