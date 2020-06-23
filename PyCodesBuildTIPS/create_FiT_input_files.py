#==========================================================
# Script to read IMERG data files and convert them to files
# which can be read by the FiT tracking algorithm. Converts
# rain rates to sort-of binary i.e. everything below a 
# low threshold is 0, everything between the first two 
# thresholds is 1, between second two thresholds is 2, 
# and above the last threshold is 3.
#
# James Russell 2019
#==========================================================

#==========================================================
# Import libraries
#==========================================================

# Import python libraries
from netCDF4 import Dataset
import glob
import time
import datetime
import driver_createinfiles as dc
from joblib import Parallel, delayed
import os

# Import namelist
import namelist_TIPS as nl

#==========================================================
# Initialize timer
#==========================================================

startnow = time.time()

#==========================================================
# Read data
#==========================================================

print("Generating file list")

# Generate a list of filenames with dates to search for
start = datetime.datetime.strptime(nl.starttime, "%Y%m%d")
end = datetime.datetime.strptime(nl.endtime, "%Y%m%d")
datearr = (start + datetime.timedelta(days=x) for x in range(0, (end-start).days))
filen = []
for dateobj in datearr:
  filen.append(nl.datadirin+nl.fileidin+dateobj.strftime("%Y%m%d")+"*")

# Reads directory and filenames
n = 0
for f in filen:
  if n==0:
    filenames = glob.glob(f)
  else:
    filenames = filenames+glob.glob(f)
  n=n+1

#==========================================================
# Write information
#========================================================== 

# Write filenames to a text file for reading during 
#  parallel loop
print("Writing filenames to text file")
ffn = open("filenames_ci.txt","w")
for i in range(len(filenames)):
  if i<len(filenames)-1:
    ffn.write(filenames[i]+",")
  else:
    ffn.write(filenames[i])
ffn.close()

#==========================================================
# Parallel loop over object
#========================================================== 

# Begins loop
if nl.serialorparallelc==1:
  print("Begin serial loop over objects")
  for i in range(len(filenames)):
    dc.driver_createinfiles(i)

# Parallel loop over PFs
if nl.serialorparallelc==2:
  print("Begin parallel loop over objects")
  Parallel(n_jobs=nl.njobsc)(delayed(
    dc.driver_createinfiles)(i) for \
    i in range(len(filenames)))

#==========================================================
# Final clean up tasks
#==========================================================

# Remove namelist and filename files
print("Cleaning up")
os.system("rm filenames_ci.txt")
os.system("rm -rf __pycache__")

# End timer
endnow = time.time()
print("Program took "+str(endnow - startnow)+"s")

#==========================================================
# End code
#==========================================================

