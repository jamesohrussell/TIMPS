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
import driver_addvars as da
from joblib import Parallel, delayed
import TIPS_functions as fns
import time as tm
import os

#==================================================================
# Namelist
#==================================================================

# Directory and filename for PF files
datadir = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_013120/IMERGPFs_2018/"
fileid  = "IPFS_"

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

# Variables desired
addmaxrr          = False # Maximum rain rate
addmeanrr         = False # Mean rain rate
addmedianrr       = False # Median rain rate
addstddevrr       = False # Standard deviation of the rain
                          #  rates
addarea           = False # Area of the PF
addvrr            = False # Volumetric rain rate
addpropagation    = False # Propagation characteristics
addnormtime       = False # Add a normalized time variable 
                          #  (0-1 for PF)
addTCinfo         = False # Flags indicating proximity to 
                          #  TC center
addlandinfo       = False # Flags indicating locations over 
                          #  land
addboundaryinfo   = False # Time-series indicating if PF 
                          #  touching domain boundary
addlocaltime      = False # Local solar time of the PF
addasymmetry      = False # Asymmetry shape parameter 
                          #  (Zick et al. 2016)
addfragmentation  = False # Fragmentation shape parameter 
                          #  (Zick et al. 2016)
addaxesshape      = False # Array of variables based on 
                          #  the major and minor axes from 
                          #  eigenvalue/vectors
addperimeter      = False # Distance around perimeter of 
                          #  shape (alpha-shape method)
addconvrain       = False # Flags indicating whether rain 
                          #  is convective
addconvarea       = False # Area of the convective region 
                          #  (addconvrain must also be True)
addconvvrr        = False # Volume of convective rainfall 
                          #  (addconvrain must also be True)
# Thermodynamic variables
addCAPEE5         = False # ERA5 Convective Available 
                          #  Potential Energy

# Moisture variables
addTCWVE5         = False # ERA5 Total Column Water Vapor
addSPHFE5         = True  # ERA5 850-200 hPa specific humidity
                          #  (free troposphere)
addSPHBE5         = False # ERA5 1000-850 hPa specific humidity 
                          #  (boundary layer)

# Kinematic/dynamics variables
addSHRFE5         = False # ERA5 850-200 hPa wind shear 
                          #  (free troposphere)
addSHRBE5         = False # ERA5 1000-850 hPa wind shear 
                          #  (boundary layer)

# Cloud variables
addCCTOE5         = False # ERA5 Total Cloud Cover


#addconvshape      = False # All shape parameters specified 
                           #  but for convective region
#addmergeinfo      = False # Information on mergers

# Inputs for specific variables

# Grid spacing in degrees lon, lat (only required for 
#  area/volrainrate/shape)
dx = 0.1
dy = 0.1

# Directory and filename of TC data (only required for TC 
#  information)
dataTCdir = "/uufs/chpc.utah.edu/common/home/varble-group2/IBTrACS/"
fileTCid  = "IBTrACS.ALL.v04r00.nc"

# Convective rain rate threshold (mm/hr)
convrainthold = 10

# Directory and filenames of ERA5 data
dataE5dir    = "/uufs/chpc.utah.edu/common/home/varble-group1/ERA5/"
hda          = 5 # Half data area in degrees
fileCAPEE5id = "ERA5.CAPE."
fileTCWVE5id = "ERA5.TCWV."
fileCCTOE5id = "ERA5.CCTO."
fileSHRFE5id = "ERA5.WSHR_850-200hPamean."
fileSHRBE5id = "ERA5.WSHR_1000-850hPamean."
fileSPHFE5id = "ERA5.SPHU_850-200hPamean."
fileSPHBE5id = "ERA5.SPHB_1000-850hPamean."

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
namelist["addmaxrr"] = str(addmaxrr)
namelist["addmeanrr"] = str(addmeanrr)
namelist["addmedianrr"] = str(addmedianrr)
namelist["addstddevrr"] = str(addstddevrr)
namelist["addarea"] = str(addarea)
namelist["addvrr"] = str(addvrr)
namelist["addpropagation"] = str(addpropagation)
namelist["addnormtime"] = str(addnormtime)
namelist["addTCinfo"] = str(addTCinfo)
namelist["addlandinfo"] = str(addlandinfo)
namelist["addboundaryinfo"] = str(addboundaryinfo)
namelist["addlocaltime"] = str(addlocaltime)
namelist["addconvrain"] = str(addconvrain)
namelist["addconvarea"] = str(addconvarea)
namelist["addconvvrr"] = str(addconvvrr)
namelist["addaxesshape"] = str(addaxesshape)
namelist["addperimeter"] = str(addperimeter)
namelist["addasymmetry"] = str(addasymmetry)
namelist["addfragmentation"] = str(addfragmentation)
namelist["addCAPEE5"] = str(addCAPEE5)
namelist["addTCWVE5"] = str(addTCWVE5)
namelist["addCCTOE5"] = str(addCCTOE5)
namelist["addSHRBE5"] = str(addSHRBE5)
namelist["addSHRFE5"] = str(addSHRFE5)
namelist["addSPHBE5"] = str(addSPHBE5)
namelist["addSPHFE5"] = str(addSPHFE5)

if addarea or addvrr or addconvarea or addconvvrr or \
   addperimeter or addasymmetry or addfragmentation or \
   addaxesshape or addboundaryinfo:
  namelist["dx"] = dx
  namelist["dy"] = dy
if addTCinfo:
  namelist["dataTCdir"] = str(dataTCdir)
  namelist["fileTCid"] = str(fileTCid)
if addconvrain or addconvarea or addconvvrr:
  namelist["convrainthold"] = convrainthold
if addCAPEE5 or addTCWVE5 or addCCTOE5:
  namelist["dataE5dir"] = str(dataE5dir)
  namelist["hda"] = hda
if addCAPEE5: namelist["fileCAPEE5id"] = str(fileCAPEE5id)
if addTCWVE5: namelist["fileTCWVE5id"] = str(fileTCWVE5id)
if addCCTOE5: namelist["fileCCTOE5id"] = str(fileCCTOE5id)
if addSHRBE5: namelist["fileSHRBE5id"] = str(fileSHRBE5id)
if addSHRFE5: namelist["fileSHRFE5id"] = str(fileSHRFE5id)
if addSPHBE5: namelist["fileSPHBE5id"] = str(fileSPHBE5id)
if addSPHFE5: namelist["fileSPHFE5id"] = str(fileSPHFE5id)

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
