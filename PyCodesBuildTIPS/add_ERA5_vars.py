#==================================================================
# Add ERA5 variables to IMERG PF netcdf files. Requires 
#  driver_addERA5vars.py be in the same directory.
#
# James Russell 2019
#==================================================================

# Import python libraries (do not change)
import numpy as np
from netCDF4 import Dataset
import glob
import datetime as dt
import driver_addERA5vars as da
from joblib import Parallel, delayed
import time as tm
import os

# Import namelist
from namelist_TIPS import general as gnl
from namelist_TIPS import addERA5vars as anl

#==================================================================
# Initialize timer
#==================================================================

startnow = tm.time()

#==================================================================
# Generate list of IPF files to process
#==================================================================

# Reads directory and file names
if anl.ssdat:

  print("Subsetting by date")

  # Generate a list of filenames with dates to search for
  print("Generating filenames to search for")
  start = dt.datetime.strptime(anl.date1,"%Y%m%d")
  end = dt.datetime.strptime(anl.date2,"%Y%m%d")
  datearr = (start + dt.timedelta(days=x) for x in range(0,(end-start).days))
  filen = []
  for dateobj in datearr:
    filen.append(gnl.datadirTIPS+gnl.fileidTIPS+"*"+dateobj.strftime("%Y%m%d")+"*")

  # Reads directory and filenames
  print("Getting directory and filenames")
  filenamesrun = []
  for f in filen:
      filenamesrun = filenamesrun+sorted(glob.glob(f))

elif anl.ssobs:

  print("Subsetting by object ids")

  print("Generating filenames to search for")
  obids = [i for i in range(int(anl.obid1),int(anl.obid2))]
  filen = []
  for i in range(len(obids)):
    filen.append(gnl.datadirTIPS+gnl.fileidTIPS+str(obids[i])+"*")

  # Reads directory and filenames
  print("Getting directory and filenames")
  filenamesrun = []
  for f in filen:
    print("Adding "+str(f))
    filenamesrun = filenamesrun+sorted(glob.glob(f))

else:

  print("No subsetting")

  # Find all files in directory
  print("Generating list of files")
  filenamesrun = sorted(glob.glob(gnl.datadirTIPS+gnl.fileidTIPS+'*'))

# Write file names to a text file for reading during 
#  parallel loop
print("Writing filenames to text file")
ffn = open("filenames_av_E5.txt","w")
for i in range(len(filenamesrun)):
  ffn.write(filenamesrun[i]+"\n")
ffn.close()

#==================================================================
# Loop over IPF files in parrallel
#==================================================================

if anl.serialorparallel==1:
  print("Begin serial loop over TIPS")
  for fn in range(len(filenamesrun)):
    da.driver_addvars(fn)

# Parrallel loop over PFs
if anl.serialorparallel==2:
  print("Begin parrallel loop over TIPS")
  Parallel(n_jobs=anl.njobs)(delayed(da.driver_addvars)(fn) \
    for fn in range(len(filenamesrun)))

#==================================================================
# Final clean up tasks
#==================================================================

# Remove namelist and filename files
print("Cleaning up")
os.system("rm filenames_av_E5.txt")
os.system("rm -rf __pycache__")

# End timer
endnow = tm.time()
print("Program took "+str(endnow - startnow)+"s")

#==================================================================
# End code
#==================================================================
