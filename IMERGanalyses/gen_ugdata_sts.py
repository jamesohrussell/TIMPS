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
fileidin  = "fn_ug_v3_"
startdate = "20180601"
enddate   = "20181001"

# For output
datadirout = ""
fileidout  = "ugdata_sts_2018_0-20N_30W-20W.nc"

latN = 20
latS = 0
lonE = -20
lonW = -30

#==========================================================
# Get universal information
#==========================================================

# Get list of filenames 
start = dt.datetime.strptime(startdate, "%Y%m%d")
end = dt.datetime.strptime(enddate, "%Y%m%d")
filenames = []
years = np.arange(int(startdate[0:4]),int(enddate[0:4])+1,1)
for iy in years:
  filenames.extend(genfromtxt(datadirin+fileidin+str(iy)+".txt", delimiter=',',dtype=str))

# Generate grid
time1 = float(start.strftime("%j"))
time2 = float(end.strftime("%j"))
tims1 = np.arange(time1,time2,1.)
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
  volr    = datasetFin.variables["volrainrate"]
  area    = datasetFin.variables["area"]

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
    datetim1  = [i for i in datetim[cmninds]]

    # Find index of interval it falls within
    for it in range(len(centlon1)):

      if centlat1[it]>latN:
        print("Ignoring - outside domain")
        continue
      elif centlat1[it]<latS:
        print("Ignoring - outside domain")
        continue
      elif centlon1[it]>lonE:
        print("Ignoring - outside domain")
        continue
      elif centlon1[it]<lonW:
        print("Ignoring - outside domain")
        continue

    # Get current time as day of year
    timenow = dt.datetime.strptime(str(datetim1[it]), "%Y%m%d%H%M")
    timey   = float(timenow.strftime("%j"))

    # Find index of interval it falls within
    idt = (np.abs(tims1-timey)).argmin()

    # Add to bin
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



