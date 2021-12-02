#=====================================================================
# Script to add variables to track files
#
# James Russell 2020
#=====================================================================

# Import python libraries (do not change)
import numpy as np
import glob
import datetime
from netCDF4 import Dataset
import os
import sys
from joblib import Parallel, delayed
import time

sys.path.insert(0,
'/uufs/chpc.utah.edu/common/home/u0816744/general_functions/')
import misc_functions as mfns

import warnings
warnings.filterwarnings("ignore")

# Import namelist
from namelist_TIMPS import general as gnl

#=====================================================================
# Initialize timer
#=====================================================================

startnow = time.time()

#=====================================================================
# Get all IDs that are MCSs
#=====================================================================

print("Getting all MCS IDs")
MCSfiles = sorted(glob.glob(f'{gnl.datadirTIMPS}/*/\
{gnl.fileidTIMPS}*'))
inlen = len(gnl.datadirTIMPS)+3+len(gnl.fileidTIMPS)
MCSIDs = np.array([int(i[inlen:inlen+7]) for i in MCSfiles])

#=====================================================================
# Get all tracking files

Trkfiles = sorted(glob.glob(f'{gnl.datadirtrkout}/*/\
{gnl.fileidtrkout}*'))

#=====================================================================
# Begin loop over each Track file
#=====================================================================

def loop(k):

  warnings.filterwarnings("ignore")

  print(f"Working on file: {Trkfiles[k]}")

#=====================================================================
# Read in track file and get IDs

  # Open file in append mode and read data
  ftr  = Dataset(Trkfiles[k])
  SID = ftr.variables["SystemID"][:]
  ftr.close()

  # Get all IDs present at this time
  SIDsTrk = np.array([int(i) for i in np.unique(SID)])

  # Get all MCS IDs present at this time
  MCSIDsnow = np.intersect1d(SIDsTrk,MCSIDs,assume_unique=True)

  # Get all non-MCS IDs present at this time
  nonMCSIDsnow = set(SIDsTrk)-set(MCSIDs)

#=====================================================================
# Build SID array with only MCSs

  # Create zero array like SID
  SIDMCS = np.zeros_like(SID)

  # If there are MCSs, add them in
  if len(MCSIDsnow)>0:
    SIDMCS = np.where(np.isin(SID, MCSIDsnow, invert=True), 0, SID)

#=====================================================================
# Write to file

  # Open file in append mode
  fta  = Dataset(Trkfiles[k],"a")

  # Create variable in file
  try:
    SystemID_MCS  = fta.createVariable('SystemID_MCS',
     np.float64,('time','lat','lon'),zlib=True,complevel=1)
  except:
    SystemID_MCS = fta.variables["SystemID_MCS"]

  # Write variable
  SystemID_MCS[:] = SIDMCS
 
  # Close file
  fta.close()

#=====================================================================
# Parallel loop over object
#===================================================================== 

if gnl.wtserialorparallel==1:
  print("Starting serial loop")
  for k in range(len(Trkfiles)):
    loop(k)

elif gnl.wtserialorparallel==2:
  print("Starting parallel loop")
  Parallel(n_jobs=gnl.wtnjobs)(delayed(loop)(k) for \
    k in range(len(Trkfiles)))

#=====================================================================
# Final clean up tasks
#=====================================================================

# Remove namelist and filename files
print("Cleaning up")
os.system("rm -rf __pycache__")

# End timer
endnow = time.time()
print("Program took "+str(endnow - startnow)+"s")

#=====================================================================
# End code
#=====================================================================







