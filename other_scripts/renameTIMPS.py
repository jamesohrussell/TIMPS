import glob
import datetime as dt
from joblib import Parallel, delayed
import os
import sys
from netCDF4 import Dataset

filenames = glob.glob("TIMPS_2015/TIMPS_*_*.nc")

# Define function
def loop(f):

  filenamenow = filenames[f]
  print(f'Working on {filenamenow}')

  # Get fn system ID
  fnID = str(filenamenow[17:24]).zfill(7)

  # Get actual system ID  
  filenow = Dataset(filenamenow,"a")
  filenow.renameAttribute("objectID","SystemID")
  filenow.close()
#  filenow = Dataset(filenamenow)
#  actualID = str(filenow.SystemID).zfill(7)
#  filenow.close()

#  # Check if they are equal
#  if not fnID==actualID:
#    newname = f'{filenamenow[0:17]}{actualID}{filenamenow[24:]}'
#    print(f'mv {filenamenow} {newname}')
#    os.system(f'mv {filenamenow} {newname}')
#    exit()

# Parallel loop over PFs
print("Begin parallel loop over objects")
Parallel(n_jobs=int(20))(delayed(
 loop)(f) for f in range(len(filenames)))
