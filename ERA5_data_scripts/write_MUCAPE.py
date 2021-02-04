# Calculate and write most unstable CAPE for ERA5

# Import libraries
import xarray as xr
import numpy as np
from netCDF4 import Dataset
import sys
sys.path.insert(0,'/uufs/chpc.utah.edu/common/home/u0816744/general_functions/')
import misc_functions as mfns
import time_functions as tfns
from joblib import Parallel, delayed
import time as tm
import glob
import sharppy.sharptab.profile as profile

# Namelist
datadirout = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/CAPE/"
dirTin     = "/uufs/chpc.utah.edu/common/home/u0816744/ERA5/temperature/allplevs/"
fileidTin  = "ERA5.TEMP."
dirQin     = "/uufs/chpc.utah.edu/common/home/u0816744/ERA5/moisture/SPHU/allplevs/"
fileidQin  = "ERA5.SPHU."
dirGin     = "/uufs/chpc.utah.edu/common/home/u0816744/ERA5/geopotential/allplevs/"
fileidGin  = "ERA5.GEOP."
dirGsin  = "/uufs/chpc.utah.edu/common/home/u0816744/ERA5/geopotential/sfc/"
fileidGsin = "ERA5.GEOP_sfc."
dirPsin  = "/uufs/chpc.utah.edu/common/home/u0816744/ERA5/pressure/sfc/"
fileidPsin = "ERA5.PRES_sfc."

# Find all temperature files
print("Finding all files to process")
Tfiles = sorted(glob.glob(dirTin+fileidTin+"*"))

def write_calc(n):

  print("Working on file "+Tfiles[n])

  # Open temperature
  print("Get temperature")
  T = xr.open_dataset(Tfiles[n])
  try:
    T = T.drop("utc_date").to_array().squeeze(
      ).sel(latitude=slice(35,-35))
  except:
    T = T.to_array().squeeze(
      ).sel(latitude=slice(35,-35))

  # Get coordinates
  time = T.time.values
  tnum = [tfns.time_since(t[0:4]+"-"+t[5:7]+"-"+t[8:10]+" "+\
                          t[11:13]+":"+t[14:16]+":"+t[17:19],
           "hours since 1900-01-01 00:00:00")
           for t in np.datetime_as_string(time)]
  levs = T.level.values
  lats = T.latitude.values
  lons = T.longitude.values

  # Generate data array for pressure
  print("Get pressure")
  P_int = T.coords['level'].values
  P = np.full(T.shape, P_int.reshape(1, -1, 1, 1))
  P = xr.DataArray(P, dims = ("time","level","latitude","longitude"),
    coords={"time":time,"level":levs,
            "latitude":lats,"longitude":lons})

  # Read specific humidity and calculate mixing ratio
  print("Get mixing ratio")
  Q = xr.open_dataset(dirQin+fileidQin+\
   Tfiles[n][len(dirTin)+len(fileidTin):-3]+".nc")
  try:
    Q = Q.drop("utc_date").sel(latitude=slice(35,-35))
  except:
    Q = Q.sel(latitude=slice(35,-35))
  w = Q/(1-Q)
  w = w.rename({'Q':'w'}).to_array().squeeze()

  # Read geopotential
  print("Get geopotential")
  G = xr.open_dataset(dirGin+fileidGin+\
   Tfiles[n][len(dirTin)+len(fileidTin):-3]+".nc")
  try:
    G = G.drop("utc_date").to_array().sel(
      latitude=slice(35,-35)).squeeze()/9.80665
  except:
    G = G.to_array().sel(
      latitude=slice(35,-35)).squeeze()/9.80665

  # Read surface pressure
  print("Get surface pressure")
  ysearch = np.datetime_as_string(time[0])[0:4]
  msearch = np.datetime_as_string(time[0])[5:7]
  dsearch = np.datetime_as_string(time[0])[8:10]
  filePs = glob.glob(dirPsin+fileidPsin+str(ysearch)+str(msearch)+"*")
  Ps = xr.open_dataset(filePs[0])
  try:
    Ps = Ps.drop("utc_date").to_array().sel(
      time=slice(ysearch+"-"+msearch+"-"+dsearch+" 00:00:00",
                 ysearch+"-"+msearch+"-"+dsearch+" 23:00:00"),
      latitude=slice(35,-35)).squeeze()/100.
  except:
    Ps = Ps.to_array().sel(
      time=slice(ysearch+"-"+msearch+"-"+dsearch+" 00:00:00",
                 ysearch+"-"+msearch+"-"+dsearch+" 23:00:00"),
      latitude=slice(35,-35)).squeeze()/100.

  # Read surface geopotential and calculate terrain height
  print("Get surface height")
  Zt = xr.open_dataset(dirGsin+fileidGsin+"nc")
  try:
    Zt = Zt.drop("utc_date").to_array().sel(
      latitude=slice(35,-35)).squeeze()
  except:
    Zt = Zt.to_array().sel(
      latitude=slice(35,-35)).squeeze()
  Zt = Zt/9.80665
  Zt = Zt.values
  Zt = np.full(Ps.shape,np.expand_dims(Zt,0)) 
  Zt = xr.DataArray(Zt, dims = ("time","latitude","longitude"),
    coords={"time":time,"latitude":lats,"longitude":lons})

  # Do CAPE calculations
  print("Calculate CAPE etc.")
  conv = wrf.cape_2d(P,T,w,G,Zt,Ps,False,meta=False)
  CAPE = conv[0,:,:,:]
  CIN  = conv[1,:,:,:]
  LCL  = conv[2,:,:,:]
  LFC  = conv[3,:,:,:]

  # Creating netcdf file to write to
  dates = Tfiles[n][len(dirTin)+len(fileidTin):-3]
  print("Writing data for time period :"+dates)
  filename = "ERA5.MUCAPE."+dates+".nc"

  print("Creating file: "+datadirout+filename)
  fileout = Dataset(datadirout+filename,"w",format="NETCDF4")

  # Global Attributes
  fileout.long_name   = "Most unstable CAPE and associated parameters"
  fileout.description = ""
  fileout.history     = 'Created ' + tm.ctime(tm.time())
  fileout.source      = 'ERA5: NCAR RDA'

  # Create dimensions in file
  time = fileout.createDimension('time',None)
  lat  = fileout.createDimension('latitude',len(lats))
  lon  = fileout.createDimension('longitude',len(lons))

  # Create dimension variables in file
  mfns.write_var("time","Time","","time",np.float64,
                 "hours since 1900-01-01 00:00:00",
                 fileout,tnum)

  mfns.write_var("latitude","Latitude","","latitude",np.float64,
                "degreesN",fileout,lats)

  mfns.write_var("longitude","Longitude","","longitude",np.float64,
                "degreesE",fileout,lons)

  # Create variables in file
  mfns.write_var("CAPE",
   "Most Unstable Convective Available Potential Energy","",
   ("time","latitude","longitude"),np.float64,"J kg**-1",fileout,CAPE)

  # Create variables in file
  mfns.write_var("CIN","Convective inhibition",
   "CIN associated with most unstable CAPE",
   ("time","latitude","longitude"),np.float64,"J kg**-1",fileout,CIN)

  # Create variables in file
  mfns.write_var("LCL","Lifting Condensation Level",
   "LCL associated with most unstable CAPE",
   ("time","latitude","longitude"),np.float64,"m",fileout,LCL)

  # Create variables in file
  mfns.write_var("LFC","Level of Free Convection",
   "LFC associated with most unstable CAPE",
   ("time","latitude","longitude"),np.float64,"m",fileout,LFC)

  fileout.close()

# don't run more than 5 processes in parallel

# Loop over all desired months
#for n in range(len(Tfiles)):
#  write_calc(n)
#  exit()

Parallel(n_jobs=6)(delayed(write_calc)(n) 
  for n in range(len(Tfiles)))

