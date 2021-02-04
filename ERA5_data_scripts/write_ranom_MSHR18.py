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
datadir = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/shear/"
fileidinU = "ERA5.USHR_1000-850hPa."
fileidinV = "ERA5.VSHR_1000-850hPa."

datadirout = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/shear/"
fileidout1 = "ERA5.MSHR_1000-850hPa.anomaly_"

dates = ["2014-12","2020-01"] # First and last month to calculate
N = 30*24 # Running average window (30 days * 24 times in day))

# Generate months
print("Getting data list")
datelist = pd.date_range(
 dates[0],dates[1],freq='MS').strftime("%Y%m").tolist()
prevMonth = (dt.datetime.strptime(dates[0],"%Y-%m") -\
 relativedelta(months=1)).strftime(r"%Y%m")
nextMonth = (dt.datetime.strptime(dates[-1],"%Y-%m") +\
 relativedelta(months=1)).strftime(r"%Y%m")
alldl = [prevMonth] + datelist + [nextMonth]

def anomaly(d):

  print(str(d)+": Working on month: "+datelist[d])

  # Find files
  dirfilesU = [glob.glob(datadir+fileidinU+alldl[d]+"*")[0],
              glob.glob(datadir+fileidinU+alldl[d+1]+"*")[0],
              glob.glob(datadir+fileidinU+alldl[d+2]+"*")[0]]
  dirfilesV = [glob.glob(datadir+fileidinV+alldl[d]+"*")[0],
              glob.glob(datadir+fileidinV+alldl[d+1]+"*")[0],
              glob.glob(datadir+fileidinV+alldl[d+2]+"*")[0]]

  # Read data
  print(str(d)+": Reading files")
  fU0 = nc.Dataset(dirfilesU[0])
  la = fU0.variables["latitude"][:]
  lo = fU0.variables["longitude"][:]

  varU0 = fU0.variables["USHR"][:,:,:]
  fU1 = nc.Dataset(dirfilesU[1])
  varU1 = fU1.variables["USHR"][:,:,:]
  ti = fU1.variables["time"][:]
  fU2 = nc.Dataset(dirfilesU[2])
  varU2 = fU2.variables["USHR"][:,:,:]

  fV0 = nc.Dataset(dirfilesV[0])
  varV0 = fV0.variables["VSHR"][:,:,:]
  fV1 = nc.Dataset(dirfilesV[1])
  varV1 = fV1.variables["VSHR"][:,:,:]
  fV2 = nc.Dataset(dirfilesV[2])
  varV2 = fV2.variables["VSHR"][:,:,:]

  fU0.close()
  fU1.close()
  fU2.close()
  fV0.close()
  fV1.close()
  fV2.close()

  # Get shape
  shp0 = np.shape(varU0)
  shp1 = np.shape(varU1)
  shp2 = np.shape(varU2)

  # Concatenate all data in single time-series
  print(str(d)+": Concatenate files")
  varU = np.concatenate((varU0,varU1,varU2),axis=0)
  varV = np.concatenate((varV0,varV1,varV2),axis=0)

  # Calculate shear magnitude
  var = np.sqrt(varU**2 + varV**2)
  var1 = np.sqrt(varU1**2 + varV1**2)

  # Calculate rolling standard deviation
  print(str(d)+": Rolling window")
  indt = shp0[-1]
  varstd = np.zeros((len(ti),len(la),len(lo)))
  for t in range(shp1[0]):
    print(str(d)+": "+str(t))
    window = var[int(shp0[0]+t-N/2):int(shp0[0]+t+N/2)+1,:,:]
    varstd[t,:,:] = np.std(window,axis=0)
  
  # Calculate running average over window size N
  varavg = uniform_filter1d(var,size=N,axis=0)[
   shp0[0]+1:shp0[0]+shp1[0]+1,:,:]

  # Calculate anomaly
  varanom = np.array(var1 - varavg)
  varanompct = np.array((var1 - varavg)/varstd)
  
  # Creating netcdf file to write to
  filename1 = datadirout+fileidout1+datelist[d]+".nc"
  print(str(d)+": Creating file: "+filename1)
  fileout1 = nc.Dataset(filename1,"w",format="NETCDF4")

  # Global Attributes
  fileout1.history     = 'Created ' + tm.ctime(tm.time())
  fileout1.source      = 'ERA5: NCAR RDA'
  fileout1.description = 'Anomalies and averages'

  # Create dimensions in file
  time = fileout1.createDimension('time',len(ti))
  lat  = fileout1.createDimension('latitude',len(la))
  lon  = fileout1.createDimension('longitude',len(lo))

  # Create dimension variables in file
  mfns.write_var("time","Time","","time",
   np.float64,"hours since 1900-01-01 00:00:00",
   fileout1,ti)

  mfns.write_var("latitude","Latitude",
    "","latitude",np.float64,
    "degreesN",fileout1,la)

  mfns.write_var("longitude","Longitude",
    "","longitude",np.float64,
    "degreesE",fileout1,lo)

  # Create variables in file
  mfns.write_var_compress("MSHR",
   "Wind shear magnitude average over 1000-850 hPa","Full field",
   ("time","latitude","longitude"),np.float64,"kg kg**-1",
   fileout1,var1,lsd=3)

  #mfns.write_var_compress("TCWV_avg",
  # "Total column water vapor running-mean",
  # "8-week running mean",("time","latitude","longitude"),
  # np.float64,"kg kg**-1",fileout1,varavg,
  # cl=10,lsd=2)

  mfns.write_var_compress("MSHR_anom",
   "Wind shear magnitude anomaly average over 1000-850 hPa",
   "Anomaly from 30 day running mean",
   ("time","latitude","longitude"),
   np.float64,"kg kg**-1",fileout1,varanom,lsd=3)

  mfns.write_var_compress("MSHR_anomn",
   "Normalized wind shear magnitude anomaly average over 1000-850 hPa",
   "Anomaly from 30 day running mean normalized by 30 day running standard deviation ([var - var_mean]/var_std)",
   ("time","latitude","longitude"),
   np.float64,"",fileout1,varanompct,lsd=3)

  fileout1.close()
#  exit()
# Loop over all desired months
#for d in range(len(datelist)):
#  anomaly(d)

Parallel(n_jobs=4)(delayed(anomaly)(d) 
  for d in range(len(datelist)))
