#==================================================================
# Add variables to IMERG PF netcdf files. Requires 
#  driver_py be in same directory.
#
# James Russell 2019
#==================================================================

# Import python libraries (do not change)
from netCDF4 import Dataset
import glob
import datetime as dt
import time as tm
from joblib import Parallel, delayed
import os
import numpy as np
import pandas as pd
import sys
import csv

# Directory for custom functions
fnsdir = "/uufs/chpc.utah.edu/common/home/u0816744/general_functions/"

# Import custom libraries
sys.path.insert(0,fnsdir)
import time_functions as tfns
import misc_functions as mfns

# Directory and filenames
datadirTIPS   = \
"/uufs/chpc.utah.edu/common/home/zipser-group2/TIMPS/TIMPS_2015/"
fileidTIPS = "TIMPS_"

# Subset (for certain range of dates) 
ssdat = False
date1 = "20150101"
date2 = "20160101"
ssobs = False
obid1 = "1000000"
obid2 = "5000000"

#==================================================================
# Initialize timer
#==================================================================

startnow = tm.time()

#==================================================================
# Generate list of IPF files to process
#==================================================================

# Reads directory and file names
if ssdat:

  print("Subsetting by date")

  # Generate a list of filenames with dates to search for
  print("Generating filenames to search for")
  start = dt.datetime.strptime(date1,"%Y%m%d")
  end = dt.datetime.strptime(date2,"%Y%m%d")
  datearr = (start + dt.timedelta(days=x) \
   for x in range(0,(end-start).days))
  filen = []
  for dateobj in datearr:
    filen.append(f'{datadirTIPS}{fileidTIPS}*\
{dateobj.strftime("%Y%m%d")}*')

  # Reads directory and filenames
  print("Getting directory and filenames")
  filenamesrun = []
  for f in filen:
      filenamesrun = filenamesrun+sorted(glob.glob(f))

elif ssobs:

  print("Subsetting by object ids")

  print("Generating filenames to search for")
  obids = [i for i in range(int(obid1),int(obid2))]
  filen = []
  for i in range(len(obids)): filen.append(
    f'{datadirTIPS}{fileidTIPS}{str(obids[i])}*')

  # Reads directory and filenames
  print("Getting directory and filenames")
  filenamesrun = []
  for f in filen:
    print(f'Adding {str(f)}')
    filenamesrun = filenamesrun+sorted(glob.glob(f))

else:

  print("No subsetting")

  # Find all files in directory
  print("Generating list of files")
  filenamesrun = sorted(glob.glob(
   f'{datadirTIPS}{fileidTIPS}*'))

filemissing = []
for f in filenamesrun:
  print(f)
  filenow = Dataset(f)
  varlist = list(filenow.variables)
  if "area" not in varlist:
    print("Missing data")
    filemissing.append([filenow])
print(filemissing)

file = open('missing.csv', 'w+', newline ='')
with file:     
    write = csv.writer(file) 
    write.writerows(filemissing)

