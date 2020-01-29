#==========================================================
# Script reads IPF data in a folder and plots all tracks on
#   a map of the tracked domain.
#
# Must provide in namelist:
#  * Directory and filename identifiers for IPF files 
#  * Various plotting options
#
# James Russell 2019
#==========================================================

# Import python libraries (do not change)
import numpy as np
from netCDF4 import Dataset
import glob
import time
import scipy.ndimage as ndimage
import datetime
import PF_functions as PFfunc

#==========================================================
# Namelist
#==========================================================

# Directory and filename for IPFs
datadirPF = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_102519/IMERGPFs_"
fileidPF  = "IPF_"
startdate = "0601"
enddate   = "1001"
year      = "2018"

# For output
fileidout  = "fn_ug_v3_"

#==========================================================
# Initialize timer
#==========================================================

startnow = time.time()

#==========================================================
# Get universal information
#==========================================================

print("Generating file list")

# Generate a list of filenames with dates to search for
filen = []
start = datetime.datetime.strptime(year+startdate, "%Y%m%d")
end = datetime.datetime.strptime(year+enddate, "%Y%m%d")
datearr = (start + datetime.timedelta(days=x) for x in range(0, (end-start).days))
for dateobj in datearr:
  filen.append(datadirPF+year+"/"+fileidPF+"*"+dateobj.strftime("%Y%m%d")+"*")
 
# Reads directory and filenames
n = 0
for f in filen:
  if n==0:
    filenames = glob.glob(f)
  else:
    filenames = filenames+glob.glob(f)
  n=n+1

#==========================================================
# Begin loop over all IPFs
#==========================================================

filenames_ug = []

# Loop over filenames
for f in filenames:
  print("Reading "+str(f))
  datasetFin = Dataset(f)

#==========================================================
# Ensure this is a valid tropical oceanic MCS
#==========================================================

  # Check if system is associated with a TC
  if datasetFin.within_TC=="True":
    print("Ignoring object "+str(f)+" because it is a hurricane")
    continue

  # Check if system is over land
  if 1 in datasetFin.variables["cPF_over_land"]:
    print("Ignoring object "+str(f)+" because it is over land")
    continue

  # Check if system is near tracking domain boundary
  if 1 in datasetFin.variables["touchesdombound"]:
    print("Ignoring object "+str(f)+" because it is on boundary")
    continue

  # Read IPF data
  centlat = datasetFin.variables["centrallat"]
  centlon = datasetFin.variables["centrallon"]
  area    = datasetFin.variables["area"]
  volr    = datasetFin.variables["volrainrate"]
  maxr    = datasetFin.variables["maxrainrate"]

  # Check if system does not last long enough (must last 6 hours - appox  PMW overpasses)
  if len(volr)<12:
    print("Ignoring object because it does not live for long enough")
    continue

  # Check if PF reaches a minimum threshold area to be defined as an MCS
  if max(area)<2000:
    print("Ignoring object because it does not reach a threshold size for an MCS")
    continue

  # Check if PF reaches a minimum threshold rain rate to be defined as an MCS
  if max(maxr)<10:
    print("Ignoring object because it does not have any convective rainfall")
    continue

  # Find indices of maxima in times series
  idtamx = (np.abs(area-max(area))).argmin()
  idtvmx = (np.abs(volr-max(volr))).argmin()

  # Check if object maximum is at start
  if idtamx==0 or idtvmx==0:
    print("Ignoring object because it's already past it's prime")
    continue

#==========================================================
# Ensure sufficient growth is occurring
#==========================================================

  # Find indices of minimum area and volumetric rain rate in time series
  idtamn = (np.abs(area-min(area[0:idtamx+1]))).argmin()
  idtvmn = (np.abs(volr-min(volr[0:idtvmx+1]))).argmin()

  # Check if area change is large enough
  if area[idtamx]/area[idtamn]<2:
    print("Ignoring object because it's area change is small")
    continue

  # Check if VRR change is large enough
  if volr[idtvmx]/volr[idtvmn]<2:
    print("Ignoring object because it's VRR change is small")
    continue

#==========================================================
# Ensure this isn't just a single merger of systems
#==========================================================
  
  # Check if max and min in time-series are one time-step apart
  if idtamx-idtamn<=1 or idtvmx-idtvmn<=1:
    print("Ignoring object because max and min are one time-step apart")
    continue

  # Calculate increase in these factors
  dvolr = np.gradient(volr,1800)
  darea = np.gradient(area,1800)

  # Check if area change is not mostly in one time-step
  if any(darea>0.9*(area[idtamx]-area[idtamn])): 
    print("Ignoring object because majority of area change occurs in one time step")
    continue

  # Check if VRR change is not mostly in one time-step
  if any(dvolr>0.9*(volr[idtvmx]-volr[idtvmn])):
    print("Ignoring object because majority of area change occurs in one time step")
    continue

#==========================================================
# If PF passes all tests, add to list
#==========================================================

  filenames_ug.append(str(f))
      
#==========================================================
# Write filelist
#==========================================================

# Write file names to a text file for reading during 
#  parallel loop
print("Writing filenames to text file")
ffn = open(fileidout+year+".txt","w")
c=0
for f in filenames_ug:
  if c<len(filenames_ug)-1:
    ffn.write(filenames_ug[c]+",")
  else:
    ffn.write(filenames_ug[c])
  c=c+1
ffn.close()

#==========================================================
# Final clean up tasks
#==========================================================

# End timer
endnow = time.time()
print("Program took "+str(endnow - startnow)+"s")

