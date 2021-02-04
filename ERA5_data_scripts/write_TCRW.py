import netCDF4 as nc
import numpy as np
import datetime as dt
import glob
from dateutil.relativedelta import relativedelta
from scipy.ndimage.filters import uniform_filter1d
import time as tm
import sys
import pandas as pd
sys.path.insert(0,'/uufs/chpc.utah.edu/common/home/u0816744/general_functions/')
import misc_functions as mfns
import xarray as xr
from joblib import Parallel, delayed

# Namelist

# Directory and file names
datadir = "/uufs/chpc.utah.edu/common/home/varble-group1/ERA5/clouds/TCRW/"
fileidin = "ERA5.TCRW."
datadirout = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/clouds/TCRW/"
fileidout1 = "ERA5.TCRW."

# Generate months
print("Getting file list")
files = glob.glob(f'{datadir}{fileidin}*')

def process(d):

  print(str(d)+": Working on file: "+files[d])

  # Read data
  print(str(d)+": Reading files")
  f0 = nc.Dataset(files[d])
  ti = f0.variables["time"][:]
  la = f0.variables["latitude"][:]
  lo = f0.variables["longitude"][:]
  ilaN = np.where(la==35)[0][0]
  ilaS = np.where(la==-35)[0][0]
  la = la[ilaN:ilaS+1]
  var0 = f0.variables["TCRW"][:,ilaN:ilaS+1,:]
  f0.close()
   
  # Creating netcdf file to write to
  filename1 = f'{datadirout}{fileidout1}{files[d][len(datadir)+10:len(datadir)+16]}.nc'
  print(f'{str(d)}: Creating file {filename1}')
  fileout1 = nc.Dataset(filename1,"w",format="NETCDF4")

  # Global Attributes
  fileout1.history     = f'Created {tm.ctime(tm.time())}'
  fileout1.source      = 'ERA5: NCAR RDA'

  # Create dimensions in file
  time = fileout1.createDimension('time',len(ti))
  lat  = fileout1.createDimension('latitude',len(la))
  lon  = fileout1.createDimension('longitude',len(lo))

  # Create dimension variables in file
  mfns.write_var("time","Time","","time",
   np.float64,"hours since 1900-01-01 00:00:00",
   fileout1,ti,filename1)

  mfns.write_var("latitude","Latitude",
    "","latitude",np.float64,
    "degreesN",fileout1,la,filename1)

  mfns.write_var("longitude","Longitude",
    "","longitude",np.float64,
    "degreesE",fileout1,lo,filename1)

  # Create variables in file
  mfns.write_var_compress("TCRW",
   "Total column rain water","Full field",
   ("time","latitude","longitude"),np.float64,"kg m**-2",
   fileout1,var0,cl=1,lsd=2)

  fileout1.close()

# Loop over all desired months
#for d in range(len(files)):
#  process(d)

Parallel(n_jobs=16)(delayed(process)(d) 
  for d in range(len(files)))
