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
import time_functions as tfns
from joblib import Parallel, delayed
import xarray as xr
from datetime import datetime

# Namelist

# Directory and file names
datadir = "/uufs/chpc.utah.edu/common/home/zipser-group2/ERA5_derived/convergence/"
fileidin = "ERA5.MFCONV."
datadirout = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/convergence/"
fileidout1 = "ERA5.CONV_300-100hPamean.anomaly_"

dates = ["2014-12","2020-01"] # First and last month to calculate
N = 30*24 # Running average window (30 days * 24 times in day)

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
  dirfiles = sorted((glob.glob(datadir+fileidin+alldl[d]+"*")+
                     glob.glob(datadir+fileidin+alldl[d+1]+"*")+
                     glob.glob(datadir+fileidin+alldl[d+2]+"*")))

  print(str(d)+": Reading files")
  mfile = xr.open_mfdataset(dirfiles)
  la = mfile.latitude.values
  lo = mfile.longitude.values
  var = mfile.CONV31.values
  tm0 = mfile.time.squeeze().sel(time=slice(
   f'{alldl[d][0:4]}-{alldl[d][4:6]}-01',
   f'{alldl[d][0:4]}-{alldl[d][4:6]}')).values

  tm1 = mfile.time.squeeze().sel(time=slice(
   f'{alldl[d+1][0:4]}-{alldl[d+1][4:6]}-01',
   f'{alldl[d+1][0:4]}-{alldl[d+1][4:6]}')).values

  times1 = [tfns.time_since(f'{str(t)[0:10]} {str(t)[11:19]}',
   "hours since 1900-01-01 00:00:00") for t in tm1]

  # Calculate rolling standard deviation
  print(str(d)+": Rolling window")
  ltm0 = len(tm0)
  ltm1 = len(tm1)
  varstd = np.zeros((ltm1,len(la),len(lo)))
  
  # Calculate running average over window size N
  varavg = uniform_filter1d(var,size=N,axis=0)[
   ltm0+1:ltm0+ltm1+1,:,:]

  # Calculate anomaly
  var1 = var[ltm0:ltm0+ltm1,:,:]
  varanom = np.array(var1 - varavg)

  # Calculate standard deviation of anomaly
  for t in range(ltm1):
    print(str(d)+": "+str(t))
    window = abs(var[int(ltm0+t-N/2):int(ltm0+t+N/2)+1,:,:])
    varstd[t,:,:] = np.std(window,axis=0)

  # Calculate normalized anomaly
  varanompct = np.array(varanom/varstd)
  
  # Creating netcdf file to write to
  filename1 = datadirout+fileidout1+datelist[d]+".nc"
  print(str(d)+": Creating file: "+filename1)
  fileout1 = nc.Dataset(filename1,"w",format="NETCDF4")

  # Global Attributes
  fileout1.history     = 'Created ' + tm.ctime(tm.time())
  fileout1.source      = 'ERA5: NCAR RDA'
  fileout1.description = 'Anomalies and averages'

  # Create dimensions in file
  time = fileout1.createDimension('time',None)
  lat  = fileout1.createDimension('latitude',len(la))
  lon  = fileout1.createDimension('longitude',len(lo))

  # Create dimension variables in file
  mfns.write_var("time","Time","","time",
   np.float64,"hours since 1900-01-01 00:00:00",fileout1,times1)

  mfns.write_var("latitude","Latitude",
    "","latitude",np.float64,"degreesN",fileout1,la)

  mfns.write_var("longitude","Longitude",
    "","longitude",np.float64,"degreesE",fileout1,lo)

  # Create variables in file
  mfns.write_var_compress("CONV",
   "Convergence average over 300-100 hPa","Full field",
   ("time","latitude","longitude"),np.float64,"kg kg**-1",
   fileout1,var1,cl=4,lsd=5)

  #mfns.write_var_compress("CONV_avg",
  # "Total column water vapor running-mean",
  # "8-week running mean",("time","latitude","longitude"),
  # np.float64,"kg kg**-1",fileout1,varavg,
  # cl=10,lsd=2)

  mfns.write_var_compress("CONV_anom",
   "Normalized convergence average over 300-100 hPa",
   "Anomaly from 30 day running mean",
   ("time","latitude","longitude"),
   np.float64,"",fileout1,varanom,cl=4,lsd=5)

  mfns.write_var_compress("CONV_anomn",
   "Normalized convergence average over 300-100 hPa",
   "Anomaly from 30 day running mean normalized by 30 day running standard deviation ([var - var_mean]/var_std)",
   ("time","latitude","longitude"),
   np.float64,"",fileout1,varanompct,cl=4,lsd=5)

  fileout1.close()

# Loop over all desired months
#for d in range(len(datelist)):
#  anomaly(d)
#  exit()

Parallel(n_jobs=13)(delayed(anomaly)(d) 
  for d in range(len(datelist)))
