#==========================================================
# Script reads IMERG data files and tracked FiT objects,
# then reworks data such that all are on the same grid,
# and plots them with IMERG data pixels and outlines
# denoting FiT object(s) (or thresholds for testing).
#
# In development:
#  * Plot tracks of IPF centroids.
#
# Must provide in namelist:
#  * Directory and filename identifiers for IMERG files 
#      and FiT output netcdf files (the latter has no 
#      effect if only plotting IMERG)
#  * Domain bounds for data region if not global (needs to
#      be exactly the same as region tracked using FiT)
#  * Gaussian smoothing information for IMERG (should be 
#      the same as for FiT input files)
#  * Various plotting options
#
# James Russell 2019
#==========================================================

# Import python libraries (do not change)
import numpy as np
import glob
import time
import datetime
from netCDF4 import Dataset
import driver_plotIMERG as dp
from joblib import Parallel, delayed
import os

# Import namelist
import namelist_plot as nl

#==========================================================
# Initialize timer
#==========================================================

startnow = time.time()

#==========================================================
# Begin loop
#==========================================================

# Generate a list of filenames with dates to search for
print("Finding all files in date range")
start = datetime.datetime.strptime(nl.starttime, "%Y%m%d")
end = datetime.datetime.strptime(nl.endtime, "%Y%m%d")
datearr = (start + datetime.timedelta(days=x) for x in range(0, (end-start).days))
filen = []
for dateobj in datearr:
  filen.append(nl.datadirIM+nl.fileidIM+dateobj.strftime("%Y%m%d")+"*")

# Reads directory and filenames
n = 0
for f in filen:
  if n==0:
    filenames = glob.glob(f)
  else:
    filenames = filenames+glob.glob(f)
  n=n+1

# Write file names to a text file for reading during 
#  parallel loop
print("Writing filenames to text file")
ffn = open("filenames_plot.txt","w")
for i in range(len(filenames)):
  ffn.write(filenames[i]+"\n")
ffn.close()

# Begin loop over all times
#for i in range(len(filenames)):
#  print("Working on "+filenames[i])
#  dp.driver_plotIMERG(i)

# Parrallel loop over PFs
print("Begin parrallel loop over IPFs")
Parallel(n_jobs=nl.njobs)(delayed(dp.driver_plotIMERG)(i) \
  for i in range(len(filenames)))

#==========================================================
# Final clean up tasks
#==========================================================

# Remove namelist and filename files
print("Cleaning up")
os.system("rm filenames_plot.txt")
os.system("rm -rf __pycache__")

# End timer
endnow = time.time()
print("Program took "+str(endnow - startnow)+"s")

#==========================================================
# End code
#==========================================================

