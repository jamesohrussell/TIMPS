#==========================================================
# Add variables to IMERG PF netcdf files.
#
# Must provide in namelist:
#  * Directory and filename identifiers for PF files
#  * Directory and filename identifiers for IBtracs data
#     (only if TC information is desired)
#  * Variables desired
#  * Grid spacing of input data pixels in degrees 
#     (IMERG = 0.1, only required if area or volumetric 
#      rain rate is desired)
#  * Subset information if desired
#
# James Russell 2019
#==========================================================

# Import python libraries (do not change)
import numpy as np
from netCDF4 import Dataset
import glob
import datetime
import driver_addvars as da
from joblib import Parallel, delayed
import PF_functions as PFfunc
import time as tm
import os

#==========================================================
# Namelist
#==========================================================

# Directory and filename for PF files
datadir = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_012220/IMERGPFs_2018/"
fileid  = "IPF_"

# Subset (for certain range of dates) 
ssdat = False
date1 = "20140601"
date2 = "20140602"
ssobs = False
obid1 = "100000"
obid2 = "100010"

# Number of processes for parrallelization
njobs = 8

# Variables desired
addmaxrainrate    = True # Maximum rain rate
addmeanrainrate   = True # Mean rain rate
addmedianrainrate = True # Median rain rate
addstddevrainrate = True # Standard deviation of the rain rates
addarea           = True # Area of the PF
addvolrainrate    = True # Volumetric rain rate
addpropagation    = True # Propagation characteristics
addnormtime       = True # Add a normalized time variable (0-1 for PF)
addTCinfo         = True # Flags indicating proximity to TC center
addlandinfo       = True # Flags indicating locations over land
addboundaryinfo   = True # Time-series indicating if PF touching domain boundary
addlocaltime      = True # Local solar time of the PF
addasymmetry      = True # Asymmetry shape parameter (Zick et al. 2016)
addfragmentation  = True # Fragmentation shape parameter (Zick et al. 2016)
addaxesshape      = True # Array of variables based on the major and minor axes from eigenvalue/vectors
addperimeter      = True # Distance around perimeter of shape (alpha-shape method)
addconvrain       = True # Flags indicating whether rain is convective
addconvarea       = True # Area of the convective region (addconvrain must also be True)
addconvvrr        = True # Volume of convective rainfall (addconvrain must also be True)

#addconvshape      = False # All shape parameters specified but for convective region
#addmergeinfo      = False # Information on mergers
#addenvinfoERA5    = False # Various variables including moisture, shear, etc. from ERA5

# Inputs for specific variables

# Grid spacing in degrees lon, lat (only required for area/volrainrate/shape)
dx = 0.1
dy = 0.1

# Directory and filename of TC data (only required for TC information)
dataTCdir = "/uufs/chpc.utah.edu/common/home/varble-group2/IBTrACS/"
fileTCid  = "IBTrACS.ALL.v04r00.nc"

# Convective rain rate threshold (mm/hr)
convrainthold = 10

#==========================================================
# Initialize timer
#==========================================================

startnow = tm.time()

#==========================================================
# Write namelist to a dictionary
#==========================================================

# Put namelist information in dictionaries
namelist = {}
namelist["addmaxrainrate"] = str(addmaxrainrate)
namelist["addmeanrainrate"] = str(addmeanrainrate)
namelist["addmedianrainrate"] = str(addmedianrainrate)
namelist["addstddevrainrate"] = str(addstddevrainrate)
namelist["addarea"] = str(addarea)
namelist["addvolrainrate"] = str(addvolrainrate)
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
namelist["dx"] = dx
namelist["dy"] = dy
namelist["dataTCdir"] = str(dataTCdir)
namelist["fileTCid"] = str(fileTCid)
namelist["convrainthold"] = convrainthold

# Write namelist dictionary to netcdf file for reading 
#  during parallel loop
print("Writing namelist to netcdf file")
nlfileout = Dataset("namelist_av.nc","w",format="NETCDF4")
for k,v in namelist.items():
  setattr(nlfileout, k,  v)
nlfileout.close()

#==========================================================
# Generate list of IPF files to process
#==========================================================

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
      filenamesrun = filenamesrun+glob.glob(f)

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
    filenamesrun = filenamesrun+glob.glob(f)

  #import os
  #for f in filen:
  #  print(os.system("ls "+str(f)))

else:

  print("No subsetting")

  # Find all files in directory
  print("Generating list of files")
  filenamesrun = glob.glob(datadir+fileid+'*')


# Write file names to a text file for reading during 
#  parallel loop
print("Writing filenames to text file")
ffn = open("filenames_av.txt","w")
for i in range(len(filenamesrun)):
  ffn.write(filenamesrun[i]+"\n")
ffn.close()

#==========================================================
# Loop over IPF files in parrallel
#==========================================================

# Parrallel loop over PFs
print("Begin parrallel loop over IPFs")
Parallel(n_jobs=njobs)(delayed(da.driver_addvars)(fn) for \
  fn in range(len(filenamesrun)))

#==========================================================
# Final clean up tasks
#==========================================================

# Remove namelist and filename files
print("Cleaning up")
os.system("rm filenames_av.txt")
os.system("rm namelist_av.nc")
os.system("rm -rf __pycache__")

# End timer
endnow = tm.time()
print("Program took "+str(endnow - startnow)+"s")

#==========================================================
# End code
#==========================================================
