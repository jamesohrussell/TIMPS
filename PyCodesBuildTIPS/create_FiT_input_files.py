#==========================================================
# Script to read IMERG data files and convert them to files
# which can be read by the FiT tracking algorithm. Converts
# rain rates to thresholded data i.e. everything below a 
# low threshold is 0, everything between the first two 
# thresholds is 1, between second two thresholds is 2, 
# and above the last threshold is 3.
#
# Addition of Mani Rajagopal's normalized thresholds 
# methods (6/22/2019). These look at each contiguous 
# object and normalize the thresholds within each relative
# to the maximum precipation. Thus thresholds are a 
# fraction of the max precip within a contiguous object.
#
# Requires driver_createinfiles.py be in same directory.
#
# James Russell 2019
#==========================================================

#==========================================================
# Import libraries
#==========================================================

# Import python libraries
import glob
import time
import datetime
from joblib import Parallel, delayed
import os

# Import namelist
from namelist_TIPS import general as gnl
from namelist_TIPS import create as cnl

# Import scripts
sys.path.insert(0,gnl.scriptsdir)
import driver_createinfiles as dc

#==========================================================
# Initialize timer
#==========================================================

startnow = time.time()

#==========================================================
# Read data
#==========================================================

print("Generating file list")

# Generate a list of filenames with dates to search for
start = datetime.datetime.strptime(gnl.starttime, "%Y%m%d")
end = datetime.datetime.strptime(gnl.endtime, "%Y%m%d")
datearr = (start + datetime.timedelta(days=x) for x in range(0, (end-start).days))
filen = []
for dateobj in datearr:
  filen.append(gnl.datadirin+gnl.fileidin+dateobj.strftime("%Y%m%d")+"*")

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
if cnl.serialorparallel==1:
  print("Begin serial loop over objects")
  for i in range(len(filenames)):
    dc.driver_createinfiles(i)

# Parallel loop over PFs
if cnl.serialorparallel==2:
  print("Begin parallel loop over objects")
  Parallel(n_jobs=cnl.njobs)(delayed(
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

