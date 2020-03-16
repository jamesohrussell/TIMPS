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
import datetime
import driver_addERA5vars as da
from joblib import Parallel, delayed
import TIPS_functions as fns
import time as tm
import os

#==================================================================
# Namelist
#==================================================================

# Directory and filename for PF files
datadir = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_CPEX-AW/TIPS_test/"
fileid  = "TIPS_"

# Subset (for certain range of dates) 
ssdat = False
date1 = "20180601"
date2 = "20180602"
ssobs = False
obid1 = "100000"
obid2 = "100010"

# Number of processes for parrallelization
serialorparallel = 1
njobs = 8

# Types of variables desired
addctarea         = False # Area centered on TIPS
addctmean         = True  # Mean of centered area
addctarearainchk  = False # As above but checked for pixels that have IMERG rainfall
addctmeanrainchk  = True  # As above but checked for pixels that have IMERG rainfall
addinarea         = False # All points in calculated inflow
addinmean         = False # Add mean of inflow region
addinarearainchk  = False # As above but checked for pixels that have IMERG rainfall
addinmeanrainchk  = False # As above but checked for pixels that have IMERG rainfall

# Variables desired
addCAPEE5    = True # ERA5 Convective Available 
                     #  Potential Energy
addTCWVE5    = True # ERA5 Total Column Water Vapor
addSPHFE5    = True # ERA5 850-200 hPa specific humidity
                     #  (free troposphere)
addSPHBE5    = True # ERA5 1000-850 hPa specific humidity 
                     #  (boundary layer)
addSHRFE5    = True # ERA5 850-200 hPa wind shear 
                     #  (free troposphere)
addSHRBE5    = True # ERA5 1000-850 hPa wind shear 
                     #  (boundary layer)

# Half data area in degrees
hda          = 5 

# Directory and filenames of ERA5 data
dataE5dir    = "/uufs/chpc.utah.edu/common/home/varble-group1/ERA5/"
fileCAPEE5id = "convparams/ERA5.CAPE."
fileTCWVE5id = "moisture/ERA5.TCWV."
fileTCRWE5id = "clouds/ERA5.TCRW."
fileUSRFE5id = "shear/ERA5.USHR_850-200hPamean."
fileUSRBE5id = "shear/ERA5.USHR_1000-850hPamean."
fileVSRFE5id = "shear/ERA5.VSHR_850-200hPamean."
fileVSRBE5id = "shear/ERA5.VSHR_1000-850hPamean."
fileSPHFE5id = "moisture/ERA5.SPHU_850-200hPamean."
fileSPHBE5id = "moisture/ERA5.SPHU_1000-850hPamean."

#==================================================================
# Initialize timer
#==================================================================

startnow = tm.time()

#==================================================================
# Write namelist to a dictionary
#==================================================================

# Put namelist information in dictionaries
namelist = {}

# Add which variables are selected
namelist["addCAPEE5"] = str(addCAPEE5)
namelist["addTCWVE5"] = str(addTCWVE5)
namelist["addTCRWE5"] = str(addTCRWE5)
namelist["addSHRBE5"] = str(addSHRBE5)
namelist["addSHRFE5"] = str(addSHRFE5)
namelist["addSPHBE5"] = str(addSPHBE5)
namelist["addSPHFE5"] = str(addSPHFE5)

namelist["addctarea"] = str(addctarea)
namelist["addctmean"] = str(addctmean)
namelist["addctarearainchk"] = str(addctarearainchk)
namelist["addctmeanrainchk"] = str(addctmeanrainchk)
namelist["addinarea"] = str(addinarea)
namelist["addinmean"] = str(addinmean)
namelist["addinarearainchk"] = str(addinarearainchk)
namelist["addinmeanrainchk"] = str(addinmeanrainchk)

namelist["dataE5dir"] = str(dataE5dir)

if addctarea or addctmean:
  namelist["hda"] = hda

if addCAPEE5: namelist["fileCAPEE5id"] = str(fileCAPEE5id)
if addTCWVE5: namelist["fileTCWVE5id"] = str(fileTCWVE5id)
if addSHRBE5: namelist["fileUSRBE5id"] = str(fileUSRBE5id)
if addSHRFE5: namelist["fileUSRFE5id"] = str(fileUSRFE5id)
if addSHRBE5: namelist["fileVSRBE5id"] = str(fileVSRBE5id)
if addSHRFE5: namelist["fileVSRFE5id"] = str(fileVSRFE5id)
if addSPHBE5: namelist["fileSPHBE5id"] = str(fileSPHBE5id)
if addSPHFE5: namelist["fileSPHFE5id"] = str(fileSPHFE5id)

if addctarearainchk or addctmeanrainchk or \
   addinarearainchk or addinmeanrainchk:  
  namelist["fileTCRWE5id"] = str(fileTCRWE5id)

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
  start = datetime.datetime.strptime(date1,"%Y%m%d")
  end = datetime.datetime.strptime(date2,"%Y%m%d")
  datearr = (start + datetime.timedelta(days=x) for x in range(0,(end-start).days))
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
