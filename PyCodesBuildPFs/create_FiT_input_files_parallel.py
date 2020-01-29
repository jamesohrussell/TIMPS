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

# Import python libraries
from netCDF4 import Dataset
import glob
import time
import datetime
import driver_createinfiles as dc
from joblib import Parallel, delayed
import os

#==========================================================
# Namelist
#==========================================================

# Directory and filename for input IMERG data
datadirin = "/uufs/chpc.utah.edu/common/home/varble-group2/IMERG/"
fileidin  = "3B-HHR.MS.MRG.3IMERG."
datahdf5  = True
datanc4   = False

# Date and time range
starttime = "20180601"
endtime   = "20180610" # Actual last day is day before

# Thresholds for rain rates (lowest to highest)
tholds  = [0.5,1.5,4.5,13.5]

# Directory and filename for output FiT input data
datadirout = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_012220/FiT_input_2018/"
fileidout  = "IMERG_FiT_tholds_"

# Subset regions (ranges have no affect if ssreg=False)
ssreg  = True
latN   = 30
latS   = -10
lonW   = -70
lonE   = 0
ssname = "CapVer"

# Smoothing
smooth   = True
gaussian = False
stddev   = 1.5  # Standard deviation for gaussian smoothing
uniform  = True
width    = 3    # Number of points to spread running average

# Number of cores for parallel code
njobs = 8

#==========================================================
# Write namelist to a dictionary
#==========================================================

# Put namelist information in dictionaries
namelist = {}
namelist["datadirin"] = str(datadirin)
namelist["fileidin"] = str(fileidin)
namelist["datahdf5"] = str(datahdf5)
namelist["datanc4"] = str(datanc4)
namelist["starttime"] = str(starttime)
namelist["endtime"] = str(endtime)
namelist["tholds"] = tholds
namelist["datadirout"] = str(datadirout)
namelist["fileidout"] = str(fileidout)
namelist["ssreg"] = str(ssreg)
namelist["latN"] = latN
namelist["latS"] = latS
namelist["lonE"] = lonE
namelist["lonW"] = lonW
namelist["ssname"] = str(ssname)
namelist["smooth"] = str(smooth)
namelist["gaussian"] = str(gaussian)
namelist["stddev"] = stddev
namelist["uniform"] = str(uniform)
namelist["width"] = width

# Write namelist dictionary to netcdf file for reading 
#  during parallel loop
print("Writing namelist to netcdf file")
nlfileout = Dataset("namelist_ci.nc","w",format="NETCDF4")
for k,v in namelist.items():
  setattr(nlfileout, k,  v)
nlfileout.close()

#==========================================================
# Initialize timer
#==========================================================

startnow = time.time()

#==========================================================
# Read data
#==========================================================

print("Generating file list")

# Generate a list of filenames with dates to search for
start = datetime.datetime.strptime(starttime, "%Y%m%d")
end = datetime.datetime.strptime(endtime, "%Y%m%d")
datearr = (start + datetime.timedelta(days=x) for x in range(0, (end-start).days))
filen = []
for dateobj in datearr:
  filen.append(datadirin+fileidin+dateobj.strftime("%Y%m%d")+"*")

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
#for i in range(len(filenames)):
#  dc.driver_createinfiles(i)

# Parallel loop over PFs
print("Begin parallel loop over objects")
Parallel(n_jobs=njobs)(delayed(dc.driver_createinfiles)(i) for \
  i in range(len(filenames)))

#==========================================================
# Final clean up tasks
#==========================================================

# Remove namelist and filename files
print("Cleaning up")
os.system("rm filenames_ci.txt")
os.system("rm namelist_ci.nc")
os.system("rm -rf __pycache__")

# End timer
endnow = time.time()
print("Program took "+str(endnow - startnow)+"s")

#==========================================================
# End code
#==========================================================

