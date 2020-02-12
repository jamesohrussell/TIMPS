
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
ssdat = False
date1 = "20180830"
date2 = "20180902"
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

  # Rewrite time variables to file
  datafile = fileout.variables["time"]
  datafile.long_name   = "Time"
  datafile.description = ""
  del(datafile)

  # Rewrite datetime variables to file
  datafile = fileout.variables["datetime"]
  datafile.long_name   = "Date and time"
  datafile.description = ""
  del(datafile)

  # Rewrite centrallocx variables to file
  datafile = fileout.variables["centrallocx"]
  datafile.long_name   = "Central x-location"
  datafile.description = "x-location of PF centroid in tracked domain."
  del(datafile)

  # Rewrite centrallocy to file
  datafile = fileout.variables["centrallocy"]
  datafile.long_name   = "Central y-location"
  datafile.description = "y-location of PF centroid in tracked domain."
  del(datafile)
  
  # Rewrite centrallat variables to file
  datafile = fileout.variables["centrallat"]
  datafile.long_name   = "Central latitude"
  datafile.description = "Latitude of PF centroid."
  del(datafile)

  # Rewrite centrallon variables to file
  datafile = fileout.variables["centrallon"]
  datafile.long_name   = "Central longitude"
  datafile.description = "Longitude of PF centroid."
  del(datafile)

  # Rewrite normalizedtime variables to file
  datafile = fileout.variables["normalizedtime"]
  datafile.long_name   = "Normalized time"
  datafile.description = "A fractional time with 0 the first time of the PF and 1 the last time of the PF."
  del(datafile)

  # Rewrite maxrainrate variables to file
  datafile = fileout.variables["maxrainrate"]
  datafile.long_name   = "Max rain rate"
  datafile.description = "Maximum rain rate within PF."
  del(datafile)

  # Rewrite meanrainrate variables to file
  datafile = fileout.variables["meanrainrate"]
  datafile.long_name   = "Mean rain rate"
  datafile.description = "Mean rain rate within PF excluding pixels with zero rain rate."
  del(datafile)

  # Rewrite medianrainrate variables to file
  datafile = fileout.variables["medianrainrate"]
  datafile.long_name   = "Median rain rate"
  datafile.description = "Median rain rate within PF excluding pixels with zero rain rate."
  del(datafile)

  # Rewrite stddevrainrate variables to file
  datafile = fileout.variables["stddevrainrate"]
  datafile.long_name   = "Standard deviation of rain rate"
  datafile.description = "Standard deviation of rain rate within PF excluding pixels with zero rain rate."
  del(datafile)

  # Rewrite area variables to file
  datafile = fileout.variables["area"]
  datafile.long_name   = "Area"
  datafile.description = "Area within PF excluding pixels with zero rain rate."
  del(datafile)

  # Rewrite volrainrate variables to file
  datafile = fileout.variables["volrainrate"]
  datafile.long_name   = "Volumetric rain rate"
  datafile.description = "Volumetric rain rate within PF excluding pixels with zero rain rate. Calculated as the sum over all pixels (with non-zero rain rate) of the area of the pixel multipled by the rain rate of the pixel."
  del(datafile)

  # Rewrite propspd variables to file
  datafile = fileout.variables["propspd"]
  datafile.long_name   = "Propagation speed"
  datafile.description = "Calculated as the geodesic distance travelled by centroid divided by time taken."
  del(datafile)

  # Rewrite propdir variables to file
  datafile = fileout.variables["propdir"]
  datafile.long_name   = "Propagation direction"
  datafile.description = "Calculated as the clockwise angle from north the centroid is moving toward."
  del(datafile)

  # Rewrite cPF_over_land variables to file
  datafile = fileout.variables["cPF_over_land"]
  datafile.long_name   = "Center of PF over land"
  datafile.description = "1 if center of PF is over land. 0 if not."
  del(datafile)

  if fileout.within_TC=="True":

    # Rewrite cPF_in_TC variables to file
    datafile = fileout.variables["cPF_in_TC"]
    datafile.long_name   = "Proximity of PF center to TC"
    datafile.description = "If PF center is within a maximum radius of TC center this is 1, if not, this is 0."
    del(datafile)

    # Rewrite dist_cPF_cTC variables to file
    datafile = fileout.variables["dist_cPF_cTC"]
    datafile.long_name   = "Distance of PF center to TC center"
    datafile.description = "Calculated as the geodesic distance of the PF center to TC center."
    del(datafile)

    # Rewrite TCrad_cPF variables to file
    datafile = fileout.variables["TCrad_cPF"]
    datafile.long_name   = "Radius of TC"
    datafile.description = "Calculated by finding the largest radius of a TC from IBtracs and doubling."
    del(datafile)

  # Close current file
  fileout.close()

#==========================================================
# Final clean up tasks
#==========================================================

# End timer
end = tm.time()
print("Program took "+str(end - start)+"s")

#==========================================================
# End code
#==========================================================
