y=2015
for m in $(seq -f "%02g" 1 12); do
  mkdir TIMPS_$y/$m/
  mv TIMPS_$y/TIMPS_*_$y$m*.nc TIMPS_$y/$m
done

