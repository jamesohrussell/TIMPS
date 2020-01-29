
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
import fiona
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep

#==========================================================
# Namelist
#==========================================================

# Directory and filename for PF files
datadir = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_100719/IMERGPFs/"
fileid  = "IPF_"

# Subset (for certain range of dates) 
ssdat = True
date1 = "20180828"
date2 = "20180903"
ssobs = False
obid1 = "118686"
obid2 = "185"

#==========================================================
# Initialize timer
#==========================================================

start = tm.time()

#==========================================================
# Generate list of IPF files to process
#==========================================================

# Reads directory and file names
print("Generating file list")
if ssdat:

  # Generate a list of filenames with dates to search for
  start = datetime.datetime.strptime(date1,"%Y%m%d")
  end = datetime.datetime.strptime(date2,"%Y%m%d")
  datearr = (start + datetime.timedelta(days=x) for x in range(0,(end-start).days))
  filen = []
  for dateobj in datearr:
    filen.append(datadir+fileid+"*"+dateobj.strftime("%Y%m%d")+"*")

  # Reads directory and filenames
  filenamesrun = []
  for f in filen:
      filenamesrun = filenamesrun+glob.glob(f)

elif ssobs:

  obids = [i for i in range(int(obid1),int(obid2))]
  filen = []
  for i in range(len(obids)):
    filen.append(datadir+fileid+str(obids[i])+"*")

  # Reads directory and filenames
  filenamesrun = []
  for f in filen:
      filenamesrun = filenamesrun+glob.glob(f)

else:
  # Find all files in directory
  filenamesrun = glob.glob(datadir+fileid+'*')

#==========================================================
# Loop over IPF files
#==========================================================

# Loop over PFs
print("Begin loop over IPFs")

for f in filenamesrun:
  print("Working with file "+str(f))

#==========================================================
# Write all data to file
#==========================================================

  fileout = Dataset(f,'a')

  fTC = []

  if fileout.within_TC=="True":
    print("TC: "+f)
    fTC = fTC.append(f)
    TCnames = fileout.TCname_cPF
    print(TCnames)

  # Close current file
  fileout.close()

#==========================================================
# Final clean up tasks
#==========================================================

print(fTC)

#==========================================================
# End code
#==========================================================
