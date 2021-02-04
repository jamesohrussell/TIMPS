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

import sys
import numpy as np
from netCDF4 import Dataset
import scipy.ndimage as ndimage
import h5py

# Import namelist
from namelist_TIMPS import general as gnl
from namelist_TIMPS import create as cnl

# Import scripts
sys.path.insert(0,gnl.fnsdir)
import shape_functions as sfns
import misc_functions as mfns
import time_functions as tfns

#==========================================================
# Initialize timer
#==========================================================

startnow = time.time()

#==========================================================
# Read files
#==========================================================

print("Generating file list")

# Generate a list of filenames with dates to search for
start = datetime.datetime.strptime(gnl.starttime,"%Y%m%d")
end = datetime.datetime.strptime(gnl.endtime,"%Y%m%d")
datearr = (start+datetime.timedelta(days=x) 
 for x in range(0,(end-start).days))
filen = []

#!! Change this line for different directory structure !!#
for dateobj in datearr:
  filen.append(f'{gnl.datadirin}{dateobj.strftime("%Y")}/\
{dateobj.strftime("%m")}/{gnl.fileidin}\
{dateobj.strftime("%Y%m%d")}*')

# Reads directory and filenames
n = 0
for f in filen:
  if n==0:
    filenames = glob.glob(f)
  else:
    filenames = filenames+glob.glob(f)
  n=n+1

if len(filenames)==0:
  print("No IMERG files found. Perhaps the IMERG \
directory or fileid is incorrect. Please check the \
namelist and try again.")
  exit()

#==========================================================
# Parallel loop over object
#========================================================== 

# Begins loop
print("Begin serial loop over objects")
for i in range(len(filenames)):

  # Open file
  datasetIM = h5py.File(filenames[i],'r')
  # Read variables
  rain = datasetIM["/Grid/precipitationCal"][:,:,:]

  print(np.amax(rain))
  print(np.mean(rain[rain>0.5]))
