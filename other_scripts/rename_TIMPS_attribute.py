import glob
import datetime as dt
from joblib import Parallel, delayed
import os
import sys
from netCDF4 import Dataset

filenames = sorted(glob.glob("/uufs/chpc.utah.edu/common/home/u0816744/zips2/TIMPS/TIMPS_2015/TIMPS_*_*.nc"))

# Define function
def loop(f):

  filenamenow = filenames[f]
  print(f'Working on {filenamenow}')

  # Get fn system ID
  fnID = str(filenamenow[17:24]).zfill(7)

  # Get actual system ID  
  filenow = Dataset(filenamenow,'a')
  try:
    filenow.renameAttribute("objectID","SystemID")
  except:
    print("SystemID already written")
  filenow.close()

# Parallel loop over PFs
print("Begin loop over objects")
#Parallel(n_jobs=int(20))(delayed(
# loop)(f) for f in range(len(filenames)))

for f in range(len(filenames)):
  loop(f)
