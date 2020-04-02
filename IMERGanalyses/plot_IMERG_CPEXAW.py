#==========================================================
# Script reads IMERG data files and tracked FiT objects,
# then reworks data such that all are on the same grid,
# and plots them with IMERG data pixels and black outlines
# denoting FiT object(s) (or thresholds for testing).
#
# In development:
#  * Plot tracks of IPF centroids.
#
# Must provide in namelist:
#  * Directory and filename identifiers for IMERG files 
#      and FiT output netcdf files (the latter has no 
#      effect if only plotting IMERG)
#  * Domain bounds for data region if not global (needs to
#      be exactly the same as region tracked using FiT)
#  * Gaussian smoothing information for IMERG (should be 
#      the same as for FiT input files)
#  * Various plotting options
#
# James Russell 2019
#==========================================================

# Import python libraries (do not change)
import numpy as np
import glob
import time
import datetime
from netCDF4 import Dataset
import driver_plotIMERG_CPEXAW as dp
from joblib import Parallel, delayed
import os

#==========================================================
# Namelist
#==========================================================

# Directory and filename for IMERG data (specify format --- only hdf5 or netcdf4 files)
# Should obtain from create_FiT_input_files.py
datadirIM = "/uufs/chpc.utah.edu/common/home/varble-group2/IMERG/"
fileidIM  = "3B-HHR.MS.MRG.3IMERG."
datahdf5  = True
datanc4   = False

# Date and time range
starttime = "20180701"
endtime   = "20180816" # Actual last day is day before

# Directory and filename for FiT output data files
datadirFi = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_102519/out_files/"
fileidFi1 = "IMERG_tracked_"
fileidFi2 = "_4Dobjects.nc"
fileidFi3 = "IMERG_tracked_4Dobject_tree.txt"

# Directory and filename for IPFs
#datadirPF = "/uufs/chpc.utah.edu/common/home/varble-group2/james/IMERGPrecipitationFeatures/"
#fileidPF  = "IPF_"

# Subset regions (if global set ssreg=False)
# Should obtain from create_FiT_input_files.py
ssreg  = True
latN   = 30
latS   = 0
lonW   = -45
lonE   = -5
ssname = "CPEXAW"

# Gaussian smoothing (No smoothing if False)
# Should obtain from create_FiT_input_files.py
smooth = False
gaussian = False
stddev   = 1.5  # Standard deviation for gaussian smoothing
uniform  = False
width    = 5    # Number of points to spread running average

# if plotobj=True, plots object outlines 
# if idobj=False, outlines all objects after obid1
# if idobj=True, outlines only object with obid1
# if longobj=True, plots nlong longest lived objects
plotobj = False
idobj   = False
obid1   = 100039
longobj = False
nlong   = 10

# If plotthr=True plots thresholds used for FiT
# This option should be False if plotobj = True as this 
#   will just overwrite the contours there. 
# Should obtain thresholds from create_FiT_input_files.py
plotthr = False
tholds  = [0.8,3.0,9.0,27.0]

# For plotting
plot  = True # Turn on or off plot (mostly for debugging)
labnx = 9    # Number of x-labels between/incl. lonW,lonE
labny = 7    # Number of y-labels between/incl. latS,latN
dx    = 0.1  # Grid-spacing in x (0.1 for IMERG)
dy    = 0.1  # Grid spacing in y (0.1 for IMERG)
latmn = 0    # Min lat for plotting
latmx = 26   # Max lat for plotting
lonmn = -45  # Min lon for plotting
lonmx = -5   # Max lon for plotting
mincl = -0.5 # Min for colorbar
maxcl = 16.5 # Max for colorbar
#tracks= True # Plot tracks?

# For parallelization
njobs = 8

#==========================================================
# Initialize timer
#==========================================================

startnow = time.time()

#==========================================================
# Write namelist to a dictionary
#==========================================================

# Put namelist information in dictionaries
namelist = {}
namelist["datadirIM"] = str(datadirIM)
namelist["fileidIM"] = str(fileidIM)
namelist["datahdf5"] = str(datahdf5)
namelist["datanc4"] = str(datanc4)
namelist["starttime"] = starttime
namelist["endtime"] = endtime
namelist["datadirFi"] = str(datadirFi)
namelist["fileidFi1"] = str(fileidFi1)
namelist["fileidFi2"] = str(fileidFi2)
namelist["fileidFi3"] = str(fileidFi3)
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
namelist["plotobj"] = str(plotobj)
namelist["idobj"] = str(idobj)
namelist["obid1"] = obid1
namelist["longobj"] = str(longobj)
namelist["nlong"] = nlong
namelist["plotthr"] = str(plotthr)
namelist["tholds"] = tholds
namelist["plot"] = str(plot)
namelist["labnx"] = labnx
namelist["labny"] = labny
namelist["dx"] = dx
namelist["dy"] = dy
namelist["latmn"] = latmn
namelist["latmx"] = latmx
namelist["lonmn"] = lonmn
namelist["lonmx"] = lonmx
namelist["mincl"] = mincl
namelist["maxcl"] = maxcl

# Write namelist dictionary to netcdf file for reading 
#  during parallel loop
print("Writing namelist to netcdf file")
nlfileout = Dataset("namelist_plot.nc","w",format="NETCDF4")
for k,v in namelist.items():
  setattr(nlfileout, k,  v)
nlfileout.close()

#==========================================================
# Begin loop
#==========================================================

# Generate a list of filenames with dates to search for
print("Finding all files in date range")
start = datetime.datetime.strptime(starttime, "%Y%m%d")
end = datetime.datetime.strptime(endtime, "%Y%m%d")
datearr = (start + datetime.timedelta(days=x) for x in range(0, (end-start).days))
filen = []
for dateobj in datearr:
  filen.append(datadirIM+fileidIM+dateobj.strftime("%Y%m%d")+"*")

# Reads directory and filenames
n = 0
for f in filen:
  if n==0:
    filenames = glob.glob(f)
  else:
    filenames = filenames+glob.glob(f)
  n=n+1

# Write file names to a text file for reading during 
#  parallel loop
print("Writing filenames to text file")
ffn = open("filenames_plot.txt","w")
for i in range(len(filenames)):
  ffn.write(filenames[i]+"\n")
ffn.close()

# Begin loop over all times
#for i in range(len(filenames)):
#  print("Working on "+filenames[i])
#  dp.driver_plotIMERG(i)

# Parrallel loop over PFs
print("Begin parrallel loop over IPFs")
Parallel(n_jobs=njobs)(delayed(dp.driver_plotIMERG)(i) for \
  i in range(len(filenames)))

#==========================================================
# Final clean up tasks
#==========================================================

# Remove namelist and filename files
print("Cleaning up")
os.system("rm filenames_plot.txt")
os.system("rm namelist_plot.nc")
os.system("rm -rf __pycache__")

# End timer
endnow = time.time()
print("Program took "+str(endnow - startnow)+"s")

#==========================================================
# End code
#==========================================================

