#=====================================================================
# Script to read FiT tracking output and IMERG files and combine 
#  corresponding output tracking and precipitation fields for 
#  tracking domain with corresponding lat and lon data.
#
# James Russell 2020
#=====================================================================

# Import python libraries (do not change)
import numpy as np
import glob
import datetime
from netCDF4 import Dataset
import os
import h5py
import pandas as pd
import sys
from joblib import Parallel, delayed
import time

sys.path.insert(0,
'/uufs/chpc.utah.edu/common/home/u0816744/general_functions/')
import misc_functions as mfns
import time_functions as tfns

import warnings
warnings.filterwarnings("ignore")

# Import namelist
from namelist_TIMPS import general as gnl
from namelist_TIMPS import create as cnl

#=====================================================================
# Initialize timer
#=====================================================================

startnow = time.time()

#=====================================================================
# Get list of files 
#=====================================================================

print("Getting all files")

# Get date range
Inputfiles = sorted(glob.glob(f'{gnl.datadirFiTin}\
{gnl.fileidFiTin}*'))

infile0    = Dataset(Inputfiles[0])
startdate  = infile0.datestart
enddate    = infile0.dateend
startdate1 = f'{str(startdate)[0:4]}-{str(startdate)[4:6]}-\
{str(startdate)[6:8]} 00:00:00'
enddate1   = f'{str(enddate)[0:4]}-{str(enddate)[4:6]}-\
{str(enddate)[6:8]} 00:00:00'
sdate = pd.to_datetime(startdate1)
edate = pd.to_datetime(enddate1)
dates = pd.date_range(sdate,edate,freq='30min')[:-1]

# Fit output
Outputfiles = sorted(glob.glob(f'{gnl.datadirFiTout}\
{gnl.fileidFiTout1}*.nc'))

# IMERG
IMERGfiles  = []
for i in range(len(dates)):
  date = dates[i]
  IMERGfiles.append(glob.glob(f'{gnl.datadirin}{date.year}/\
{str(int(date.month)).zfill(2)}/{gnl.fileidin}\
{date.year}{str(int(date.month)).zfill(2)}\
{str(int(date.day)).zfill(2)}-S{str(int(date.hour)).zfill(2)}\
{str(int(date.minute)).zfill(2)}00-*')[0])

# Get correct latitude and longitude data
latN    = infile0.latN
latS    = infile0.latS
lonE    = infile0.lonE
lonW    = infile0.lonW

IMERG0 = Dataset(IMERGfiles[0])
lat = IMERG0["/Grid/lat"][:]
lon = IMERG0["/Grid/lon"][:]

indlatN = mfns.k_closest(lat,latN,1)[0]
indlatS = mfns.k_closest(lat,latS,1)[0]
indlonE = mfns.k_closest(lon,lonE,1)[0]
indlonW = mfns.k_closest(lon,lonW,1)[0]

latnow = lat[indlatS:indlatN+1]
lonnow = lon[indlonW:indlonE+1]

# Get other data
tholdsin = infile0.tholds
sourcein = infile0.source
if hasattr(infile0, 'smthstddev'):
  sstd     = infile0.smthstddev
else:
  sstd     = None
if hasattr(infile0, 'smthwidth'):
  swid     = infile0.smthwidth
else:
  swid     = None

# Close files
infile0.close()
IMERG0.close()

#=====================================================================
# Function to write files
#=====================================================================

def loop(k):

  warnings.filterwarnings("ignore")

# Get all data

  date = dates[k]
  print(f'Working on date: {str(date)}')

  timenow = [tfns.time_since(str(date),
   "hours since 1900-01-01 00:00:00")]

  IMERGnow   = Dataset(IMERGfiles[k])
  infilenow  = Dataset(Inputfiles[k])
  outfilenow = Dataset(Outputfiles[k])

  pcp = np.transpose(IMERGnow["/Grid/precipitationCal"][
   0,indlonW:indlonE+1,indlatS:indlatN+1])
  pcp = pcp[None,:,:]
  tracks = outfilenow.variables['value'][0,:,:]
  tracks = tracks[None,:,:]

# Write all data

  # Creating netcdf file to write to
  fileandpath = f'{gnl.datadirtrkout}{gnl.fileidtrkout}{str(infilenow.time)}_{str(Inputfiles[k])[len(gnl.datadirFiTin)+len(gnl.fileidFiTin):len(gnl.datadirFiTin)+len(gnl.fileidFiTin)+5]}.nc'
  #print(f"Creating file: {fileandpath}")
  fileout = Dataset(fileandpath,'w',format='NETCDF4_CLASSIC')

  # Global Attributes
  fileout.description = 'TIMPS Precipitation and Tracking Fields'
  fileout.data_source = f'IMERG V06B Final: {sourcein}'

  # Smoothing attributes
  if sstd is not None:
    fileout.smoothing_command = \
     "scipy.ndimage.gaussian_filter(rain,sigma=int(smthstdd))"
    fileout.smoothing_standard_deviation = stdd
  if swid is not None:
    fileout.smoothing_command = \
     "scipy.ndimage.uniform_filter(rain,size=int(smthwidth))"
    fileout.smoothing_window_width = swid
    fileout.smoothing_window_width_units  = "pixels"
    fileout.smoothing_description = \
     "Smoothing of original IMERG data input to tracking algorithm"

  # Domain attributes
  fileout.tracking_start_date = startdate
  fileout.tracking_end_date = enddate
  if gnl.domperiodicx:
    fileout.domain_periodic_in_zonal_direction=1
  else:
    fileout.domain_periodic_in_zonal_direction=0
  fileout.domain_periodic_in_meridional_direction=0

  # Threshold attributes
  fileout.thresholds = tholdsin
  if cnl.tholdtype==3:
    fileout.thresholding_description = "Contiguous area"
  elif cnl.tholdtype==1:
    fileout.thresholding_description = "Fixed: Multiple rain rates"
    fileout.minimum_threshold = cnl.minthold
    fileout.minimum_threshold_units = "mm/hr"
  elif cnl.tholdtype==2:
    fileout.thresholding_description = "Normalized: Fractional values of the contiguous area maximum rain rate"
    fileout.minimum_threshold = cnl.minthold
    fileout.minimum_threshold_units = "mm/hr"

  # Splitting attributes
  fileout.splitting_distance   = int(gnl.splitdist)
  fileout.splitting_distance_units = "pixels"
  fileout.splitting_distance_description = \
   "Maximum distance objects can be separated but\
 remain as the same object after splitting"

  # Create dimensions in file
  time = fileout.createDimension('time', len(timenow))
  lat = fileout.createDimension('lat', len(latnow))
  lon = fileout.createDimension('lon', len(lonnow))

  # Create variables in file
  time = fileout.createVariable('time',np.float64, ('time'))
  lat  = fileout.createVariable('lat', np.float64, ('lat'))
  lon  = fileout.createVariable('lon',np.float64, ('lon'))
  Precipitation  = fileout.createVariable('Precipitation',
   np.float64, ('time','lat','lon'),zlib=True,complevel=1)
  SystemID  = fileout.createVariable('SystemID',
   np.float64, ('time','lat','lon'),zlib=True,complevel=1)

  # Write Variables
  time[:] = timenow
  lat[:] = latnow
  lon[:] = lonnow
  Precipitation[:] = pcp
  SystemID[:] = tracks

# Close all files
  IMERGnow.close()
  infilenow.close()
  outfilenow.close()

#=====================================================================
# End function
#===================================================================== 

#=====================================================================
# Parallel loop over object
#===================================================================== 

if gnl.wtserialorparallel==1:
  print("Starting serial loop")
  for k in range(len(IMERGfiles)):
    loop(k)

elif gnl.wtserialorparallel==2:
  print("Starting parallel loop")
  Parallel(n_jobs=gnl.wtnjobs)(delayed(loop)(k) for \
    k in range(len(IMERGfiles)))

#=====================================================================
# Final clean up tasks
#=====================================================================

# Remove namelist and filename files
print("Cleaning up")
os.system("rm -rf __pycache__")

# End timer
endnow = time.time()
print("Program took "+str(endnow - startnow)+"s")

#=====================================================================
# End code
#=====================================================================

