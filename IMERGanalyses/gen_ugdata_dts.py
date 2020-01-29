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
import datetime as dt
from numpy import genfromtxt

#==========================================================
# Namelist
#==========================================================

# Directory and filename for IPFs
datadirin = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_102519/IMERGanalyses/"
fileidin  = "fn_ug_v2_"
startdate = "0601"
enddate   = "1001"
startyear = 2014
endyear   = 2018

timestep = 0.5

#latN = 13
#latS = 3
#lonE = -23
#lonW = -45

# For output
datadirout = ""
fileidout  = "ugdata_dts_"+str(startyear)+"-"+str(endyear)+".nc"

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

print("Generating grids")

# Generate grids
tims1 = np.arange(0.+timestep/2.,24.+timestep/2.,timestep)
bins1 = np.zeros(len(tims1),dtype=int)

#==========================================================
# Read IPF files
#==========================================================

# Loop over filenames
for f in filenames:
  print("Reading "+str(f))
  datasetFin = Dataset(f)

  # Read IPF data
  datetim = datasetFin.variables["datetime"]
  centlon = datasetFin.variables["centrallon"]
  centlat = datasetFin.variables["centrallat"]
  area    = datasetFin.variables["area"]
  volr    = datasetFin.variables["volrainrate"]
  localdt = datasetFin.variables["localsuntime"]

  # Find indices of maxima in times series
  idtamx = (np.abs(area-max(area))).argmin()
  idtvmx = (np.abs(volr-max(volr))).argmin()

  # Calculate increase in these factors
  dvolr = np.gradient(volr,1800)
  darea = np.gradient(area,1800)

  # Find indices of minimum area and volumetric rain rate in time series
  idtamn = (np.abs(area-min(area[0:idtamx+1]))).argmin()
  idtvmn = (np.abs(volr-min(volr[0:idtvmx+1]))).argmin()

  # Find indices where all conditions are met
  indvgt0 = [i for i, val in enumerate(dvolr>0) if val]  
  indagt0 = [i for i, val in enumerate(darea>0) if val]
  indvltp = [i for i in indvgt0 if i<=idtvmx]
  indaltp = [i for i in indagt0 if i<=idtamx]
  indvall = [i for i in indvltp if i>=idtvmn]
  indaall = [i for i in indaltp if i>=idtamn]
  cmninds = list(set(indvall).intersection(set(indaall)))

  # Define the local hour and minute for the feature
  if len(cmninds)>0:

    centlon1  = [i for i in centlon[cmninds]]
    centlat1  = [i for i in centlat[cmninds]]
    localdts  = [str(i) for i in localdt[cmninds]]
    hour      = [float(i[8:10]) for i in localdts]
    minute    = [float(i[10:12]) for i in localdts]
    localstime= [hour[i]+(minute[i]/60.) for i in range(len(hour))]

    # Find index of interval it falls within
    for it in range(len(localstime)):

#      if centlat1[it]>latN:
#        print("Ignoring - outside domain")
#        continue
#      elif centlat1[it]<latS:
#        print("Ignoring - outside domain")
#        continue
#      elif centlon1[it]>lonE:
#        print("Ignoring - outside domain")
#        continue
#      elif centlon1[it]<lonW:
#        print("Ignoring - outside domain")
#        continue

      idt = (np.abs(tims1-localstime[it])).argmin()

      bins1[idt] = bins1[idt] + 1

print(bins1)

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
tims = fileout.createDimension('tims', len(tims1))

# Create variables in file
tims  = fileout.createVariable('tims' , np.float32, ('tims'))
bins  = fileout.createVariable('bins' , np.float32, ('tims'))

# Write Variables
tims[:] = tims1
bins[:] = bins1



