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
from joblib import Parallel, delayed

# Namelist

# Directory and file names
datadir = "/uufs/chpc.utah.edu/common/home/zipser-group2/ERA5_derived/moisture/RLHU/"
fileidin = "ERA5.RLHU_800-400hPamean."
datadirout = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/moisture/RLHU/"
fileidout1 = "ERA5.RLHU_800-400hPamean.anomaly_"

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
  dirfiles = [glob.glob(datadir+fileidin+alldl[d]+"*")[0],
              glob.glob(datadir+fileidin+alldl[d+1]+"*")[0],
              glob.glob(datadir+fileidin+alldl[d+2]+"*")[0]]

  # Read data
  print(str(d)+": Reading files")
  f0 = nc.Dataset(dirfiles[0])
  la = f0.variables["latitude"][:]
  lo = f0.variables["longitude"][:]
  ilaN = np.where(la==35)[0][0]
  ilaS = np.where(la==-35)[0][0]
  la = la[ilaN:ilaS+1]
  var0 = f0.variables["R"][:,ilaN:ilaS+1,:]
  f1 = nc.Dataset(dirfiles[1])
  var1 = f1.variables["R"][:,ilaN:ilaS+1,:]
  ti = f1.variables["time"][:]
  f2 = nc.Dataset(dirfiles[2])
  var2 = f2.variables["R"][:,ilaN:ilaS+1,:]

  f0.close()
  f1.close()
  f2.close()

  # Get shape
  shp0 = np.shape(var0)
  shp1 = np.shape(var1)
  shp2 = np.shape(var2)

  # Concatenate all data in single time-series
  print(str(d)+": Concatenate files")
  var = np.concatenate((var0,var1,var2),axis=0)

  # Calculate rolling standard deviation
  print(str(d)+": Rolling window")
  indt = shp0[-1]
  varstd = np.zeros((len(ti),len(la),len(lo)))
  for t in range(shp1[0]):
    print(str(d)+": "+str(t))
    window = abs(var[int(shp0[0]+t-N/2):int(shp0[0]+t+N/2)+1,:,:])
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
  mfns.write_var_compress("R",
   "Specific humidity average over 800-400 hPa","Full field",
   ("time","latitude","longitude"),np.float64,"kg kg**-1",
   fileout1,var1,lsd=5)

  #mfns.write_var_compress("TCWV_avg",
  # "Total column water vapor running-mean",
  # "8-week running mean",("time","latitude","longitude"),
  # np.float64,"kg kg**-1",fileout1,varavg,
  # cl=10,lsd=2)

  mfns.write_var_compress("R_anom",
   "Specific humidity anomaly average over 800-400 hPa",
   "Anomaly from 30 day running mean",
   ("time","latitude","longitude"),
   np.float64,"kg kg**-1",fileout1,varanom,lsd=5)

  mfns.write_var_compress("R_anomn",
   "Normalized specific humidity anomaly average over 800-400 hPa",
   "Anomaly from 30 day running mean normalized by 30 day running standard deviation ([var - var_mean]/var_std)",
   ("time","latitude","longitude"),
   np.float64,"",fileout1,varanompct,cl=4,lsd=5)

  fileout1.close()
#  exit()
# Loop over all desired months
#for d in range(len(datelist)):
#  anomaly(d)

Parallel(n_jobs=9)(delayed(anomaly)(d) 
  for d in range(len(datelist)))
