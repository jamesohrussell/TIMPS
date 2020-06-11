#==================================================================
# Add variables to IMERG PF netcdf files.
#
# Must provide in namelist:
#  * Directory and filename identifiers for PF files
#  * Directory and filename identifiers for IBtracs data
#     (only if TC information is desired)
#  * Variables desired
#  * Grid spacing of input data pixels in degrees 
#     (IMERG = 0.1, only required if area or volumetric 
#     rain rate is desired)
#  * Subset information if desired
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

#==================================================================
# Namelist
#==================================================================

# Directory for custom functions
fnsdir = "/uufs/chpc.utah.edu/common/home/u0816744/general_functions/"

# Directory and filename for PF files
datadir = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_CPEX-AW/TIPS_test2/"
fileid  = "TIPS_"

# Subset (for certain range of dates) 
ssdat = False
date1 = "20180601"
date2 = "20180602"
ssobs = False
obid1 = "103231"
obid2 = "103233"

# Number of processes for parrallelization
serialorparallel = 2
njobs = 20

# Type of area or mean
addctarea = True # Area centered on TIPS
addctmean = True  # Mean of centered area
addinarea = False # All points in calculated inflow
addinmean = False # Add mean of inflow region

# Rain check or no rain check
addrainchk   = True
addnorainchk = True

# Variables desired
addTCWVE5    = True # ERA5 Total Column Water Vapor
addCAPEE5    = True # ERA5 Convective Available Potential Energy
addSR18E5    = True # ERA5 Boundary Layer Shear
addSR82E5    = True # ERA5 Free-troposphere Shear
addSH18E5    = True # ERA5 Boundary Layer Specific Humidity
addSH82E5    = True # ERA5 Free-troposphere Specific Humidity

# ERA5 domain variables
hda           = 2.5 # Half data area in degrees
hoursbefore   = np.arange(48,0,-3) # Hours before (descending)
hoursafter    = hoursbefore[::-1] # Hours after (ascending)
avgmissfrac   = 0.8 # fraction of domain missing before average is not carried out

# Directory and filenames of ERA5 data
dataE5dir     = "/uufs/chpc.utah.edu/common/home/varble-group1/ERA5/"
fileTCWVE5id  = "moisture/ERA5.TCWV."
fileCAPEE5id  = "convparams/ERA5.CAPE."
fileUSR18E5id = "shear/ERA5.USHR_1000-850hPamean."
fileVSR18E5id = "shear/ERA5.VSHR_1000-850hPamean."
fileUSR82E5id = "shear/ERA5.USHR_850-200hPamean."
fileVSR82E5id = "shear/ERA5.VSHR_850-200hPamean."
fileSH18E5id  = "moisture/ERA5.SPHU_1000-850hPamean."
fileSH82E5id  = "moisture/ERA5.SPHU_850-200hPamean."
fileTCRWE5id  = "clouds/ERA5.TCRW."

#==================================================================
# Initialize timer
#==================================================================

startnow = tm.time()

#==================================================================
# Write namelist to a dictionary
#==================================================================

# Put namelist information in dictionaries
namelist = {}
namelist["fnsdir"] = str(fnsdir)

# Add which variables are selected
namelist["addctarea"] = str(addctarea)
namelist["addctmean"] = str(addctmean)
namelist["addinarea"] = str(addinarea)
namelist["addinmean"] = str(addinmean)

namelist["addrainchk"]   = str(addrainchk)
namelist["addnorainchk"] = str(addnorainchk)

namelist["addTCWVE5"] = str(addTCWVE5)
namelist["addCAPEE5"] = str(addCAPEE5)
namelist["addSR18E5"] = str(addSR18E5)
namelist["addSR82E5"] = str(addSR82E5)
namelist["addSH18E5"] = str(addSH18E5)
namelist["addSH82E5"] = str(addSH82E5)

namelist["dataE5dir"] = str(dataE5dir)

namelist["hda"] = hda
namelist["avgmissfrac"] = avgmissfrac
namelist["hoursbefore"] = hoursbefore
namelist["hoursafter"] = hoursafter

if addTCWVE5: namelist["fileTCWVE5id"] = str(fileTCWVE5id)
if addCAPEE5: namelist["fileCAPEE5id"] = str(fileCAPEE5id)
if addSR18E5: 
  namelist["fileUSR18E5id"] = str(fileUSR18E5id)
  namelist["fileVSR18E5id"] = str(fileVSR18E5id)
if addSR82E5: 
  namelist["fileUSR82E5id"] = str(fileUSR82E5id)
  namelist["fileVSR82E5id"] = str(fileVSR82E5id)
if addSH18E5: namelist["fileSH18E5id"] = str(fileSH18E5id)
if addSH82E5: namelist["fileSH82E5id"] = str(fileSH82E5id)
if addrainchk: namelist["fileTCRWE5id"] = str(fileTCRWE5id)

# Write namelist dictionary to netcdf file for reading 
#  during parallel loop
print("Writing namelist to netcdf file")
nlfileout = Dataset("namelist_av.nc","w",format="NETCDF4")
for k,v in namelist.items():
  setattr(nlfileout, k,  v)
nlfileout.close()

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
  datearr = (start + dt.timedelta(days=x) for x in range(0,(end-start).days))
  filen = []
  for dateobj in datearr:
    filen.append(datadir+fileid+"*"+dateobj.strftime("%Y%m%d")+"*")

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
  for i in range(len(obids)):
    filen.append(datadir+fileid+str(obids[i])+"*")

  # Reads directory and filenames
  print("Getting directory and filenames")
  filenamesrun = []
  for f in filen:
    print("Adding "+str(f))
    filenamesrun = filenamesrun+sorted(glob.glob(f))

  #import os
  #for f in filen:
  #  print(os.system("ls "+str(f)))

else:

  print("No subsetting")

  # Find all files in directory
  print("Generating list of files")
  filenamesrun = sorted(glob.glob(datadir+fileid+'*'))

# Write file names to a text file for reading during 
#  parallel loop
print("Writing filenames to text file")
ffn = open("filenames_av.txt","w")
for i in range(len(filenamesrun)):
  ffn.write(filenamesrun[i]+"\n")
ffn.close()

#==================================================================
# Loop over IPF files in parrallel
#==================================================================

if serialorparallel==1:
  print("Begin serial loop over objects")
  for fn in range(len(filenamesrun)):
    da.driver_addvars(fn)

# Parrallel loop over PFs
if serialorparallel==2:
  print("Begin parrallel loop over IPFs")
  Parallel(n_jobs=njobs)(delayed(da.driver_addvars)(fn) \
    for fn in range(len(filenamesrun)))

#==================================================================
# Final clean up tasks
#==================================================================

# Remove namelist and filename files
print("Cleaning up")
os.system("rm filenames_av.txt")
os.system("rm namelist_av.nc")
os.system("rm -rf __pycache__")

# End timer
endnow = tm.time()
print("Program took "+str(endnow - startnow)+"s")

#==================================================================
# End code
#==================================================================
