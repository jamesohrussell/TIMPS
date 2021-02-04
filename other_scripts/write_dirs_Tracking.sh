y=2015
for m in $(seq -f "%02g" 1 12); do
  mkdir Tracking_$y/$y$m/
  mv Tracking_$y/Tracking_$y$m*.nc Tracking_$y/$y$m
done

