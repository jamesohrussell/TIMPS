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
from numpy import genfromtxt

#==========================================================
# Namelist
#==========================================================

# Directory and filename for IPFs
datadirin = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_102519/IMERGanalyses/"
fileidin  = "fn_ug_"
startdate = "0601"
enddate   = "1001"
startyear = 2014
endyear   = 2018

# For grid
gridd = 0.5 # Grid delta

# For output
datadirout = ""
fileidout  = "ugdatabytime_map"+str(gridd)+"_"+str(startyear)+"-"+str(endyear)+".nc"

#==========================================================
# Initialize timer
#==========================================================

startnow = time.time()

#==========================================================
# Get universal information
#==========================================================

print("Getting file list")

# Read text file with upscale growth cases
filenames = []
years = np.arange(startyear,endyear+1,1)
for iy in years:
  filenames.extend(genfromtxt(datadirin+fileidin+str(iy)+".txt", delimiter=',',dtype=str))

# Generate a list of filenames with dates to search for
#for iy in years:
#  start = datetime.datetime.strptime(str(iy)+startdate, "%Y%m%d")
#  end = datetime.datetime.strptime(str(iy)+enddate, "%Y%m%d")
#  datearr = (start + datetime.timedelta(days=x) for x in range(0, (end-start).days))

# Find all files
#for dateobj in datearr:
#  indices = [i for i, s in enumerate(mylist) if dateobj in s]

print("Generating grids")

# Read domain
datasetFin = Dataset(filenames[0])
latN     = int(datasetFin.track_dom_latN)
latS     = int(datasetFin.track_dom_latS)
lonW     = int(datasetFin.track_dom_lonW)
lonE     = int(datasetFin.track_dom_lonE)

# Generate grid with points at center of pixels
lats1 = np.arange(latS+gridd/2,latN+gridd/2,gridd)
lons1 = np.arange(lonW+gridd/2,lonE+gridd/2,gridd)
bins1 = np.zeros((len(lats1),len(lons1)))

#==========================================================
# Check upscale growth conditions for each IPF
#==========================================================

# Loop over filenames
for f in filenames:
  print("Reading "+str(f))
  datasetFin = Dataset(f)

  # Read IPF data
  centlat = datasetFin.variables["centrallat"]
  centlon = datasetFin.variables["centrallon"]
  area    = datasetFin.variables["area"]
  volr    = datasetFin.variables["volrainrate"]

  # Find indices of maxima in times series
  idtamx = (np.abs(area-max(area))).argmin()
  idtvmx = (np.abs(volr-max(volr))).argmin()

  # Find indices of minimum area and volumetric rain rate in time series
  idtamn = (np.abs(area-min(area[0:idtamx+1]))).argmin()
  idtvmn = (np.abs(volr-min(volr[0:idtvmx+1]))).argmin()

  # Calculate increase in these factors
  dvolr = np.gradient(volr,1800)
  darea = np.gradient(area,1800)
   
  # Find indices where all conditions are met
  indvgt0 = [i for i, val in enumerate(dvolr>0) if val]  
  indagt0 = [i for i, val in enumerate(darea>0) if val]
  indvltp = [i for i in indvgt0 if i<=idtvmx]
  indaltp = [i for i in indagt0 if i<=idtamx]
  indvall = [i for i in indvltp if i>=idtvmn]
  indaall = [i for i in indaltp if i>=idtamn]
  cmninds = list(set(indvall).intersection(set(indaall)))

  # Add count to bins
  for idt in cmninds:
 
    # Find index of box it falls within
    idx   = (np.abs(lons1-centlon[idt])).argmin()
    idy   = (np.abs(lats1-centlat[idt])).argmin()

    # Add to grid bins
    bins1[idy,idx] = bins1[idy,idx] + 1#(dvolr[idt]*1800)/(volr[idtamx]-volr[idtamn])
     
#==========================================================
# Write data to netcdf file
#==========================================================

# Creating netcdf file to write to
print("Creating file: "+datadirout+fileidout)
fileout = Dataset(datadirout+fileidout,'w',format='NETCDF4_CLASSIC')

# Global Attributes
fileout.description = 'Data for plotting'
fileout.history     = 'Created ' + time.ctime(time.time())

# Create dimensions in file
lats = fileout.createDimension('lats', len(lats1))
lons = fileout.createDimension('lons', len(lons1))

# Create variables in file
lats  = fileout.createVariable('lats' , np.float32, ('lats'))
lons  = fileout.createVariable('lons' , np.float32, ('lons'))
bins  = fileout.createVariable('bins' , np.float32, ('lats','lons'))

# Write Variables
lats[:] = lats1
lons[:] = lons1
bins[:,:] = bins1

#==========================================================
# Final clean up tasks
#==========================================================

# End timer
endnow = time.time()
print("Program took "+str(endnow - startnow)+"s")

