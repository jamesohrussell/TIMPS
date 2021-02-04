# Calculate and write most unstable CAPE for ERA5

# Import libraries
import xarray as xr
import numpy as np
import wrf
from netCDF4 import Dataset
import sys
sys.path.insert(0,'/uufs/chpc.utah.edu/common/home/u0816744/general_functions/')
import geophys_functions as gfns
import misc_functions as mfns
import time_functions as tfns
from joblib import Parallel, delayed
import time as tm
import glob

# Namelist
datadirout = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/convergence/"
dirQin     = "/uufs/chpc.utah.edu/common/home/u0816744/ERA5/moisture/SPHU/allplevs/"
fileidQin  = "ERA5.SPHU."
dirUin = "/uufs/chpc.utah.edu/common/home/varble-group1/ERA5/winds/UWND/allplevs/"
dirVin = "/uufs/chpc.utah.edu/common/home/varble-group1/ERA5/winds/VWND/allplevs/"
fileidUin  = "ERA5.UWND."
fileidVin  = "ERA5.VWND."

# Find all temperature files
print("Finding all files to process")
files = sorted(glob.glob(dirQin+fileidQin+"*"))

def write_calc(n):

  print("Working on file "+files[n])

  # Open temperature
  print("Get humidity")
  Q = xr.open_dataset(files[n])
  try:
    Q = Q.drop("utc_date").to_array().squeeze(
      ).sel(latitude=slice(35,-35),level=slice(100,1000))
  except:
    Q = Q.to_array().squeeze(
      ).sel(latitude=slice(35,-35),level=slice(100,1000))

  # Get coordinates
  time = Q.time.values
  tnum = [tfns.time_since(t[0:4]+"-"+t[5:7]+"-"+t[8:10]+" "+\
                          t[11:13]+":"+t[14:16]+":"+t[17:19],
           "hours since 1900-01-01 00:00:00")
           for t in np.datetime_as_string(time)]
  levs = Q.level.values
  lats = Q.latitude.values
  lons = Q.longitude.values
  dx = [gfns.calc_distance(lons[0],l,lons[1],l) for l in lats]
  dy = gfns.calc_distance(lons[0],lats[0],lons[0],lats[1])

  # Read geopotential
  print("Get U")
  U = xr.open_dataset(dirUin+fileidUin+\
   files[n][len(dirQin)+len(fileidQin):-3]+".nc")
  try:
    U = U.drop("utc_date").to_array().sel(
      latitude=slice(35,-35),level=slice(100,1000)).squeeze()
  except:
    U = U.to_array().sel(
      latitude=slice(35,-35),level=slice(100,1000)).squeeze()

  # Read geopotential
  print("Get V")
  V = xr.open_dataset(dirVin+fileidVin+\
   files[n][len(dirQin)+len(fileidQin):-3]+".nc")
  try:
    V = V.drop("utc_date").to_array().sel(
      latitude=slice(35,-35),level=slice(100,1000)).squeeze()
  except:
    V = V.to_array().sel(
      latitude=slice(35,-35),level=slice(100,1000)).squeeze()
    
  # Calculate moisture convergence
  print("Calculating divergence")
  shp = np.shape(U)
  dudx = np.transpose(np.array([np.gradient(U[:,:,i,:],dx[i],axis=2)
   for i in range(shp[2])]),(1,2,0,3))
  dvdy = np.gradient(V,dy,axis=2)
  conv  = -dudx + dvdy
  Qconv = Q*conv

  # Calculate moisture advection
  print("Calculating advection")
  dQdx = np.transpose(np.array([np.gradient(Q[:,:,i,:],dx[i],axis=2)
   for i in range(shp[2])]),(1,2,0,3))
  dQdy = np.gradient(Q,dy,axis=2)
  Qadv  = -U*dQdx + V*dQdy

  # Calculate moisture flux divergence
  print("Calculating mfc")
  Qfc = Qadv+Qconv

  # Do vertical averaging
  conv18 = np.mean(conv[:,20:,:,:],axis=1)
  conv31 = np.mean(conv[:,:8,:,:],axis=1)
  Qadv18 = np.mean(Qadv[:,20:,:,:],axis=1)
  Qadv84 = np.mean(Qadv[:,9:19,:,:],axis=1)
  Qconv18 = np.mean(Qconv[:,20:,:,:],axis=1)
  Qconv84 = np.mean(Qconv[:,9:19,:,:],axis=1)
  Qfc18 = np.mean(Qfc[:,20:,:,:],axis=1)
  Qfc84 = np.mean(Qfc[:,9:19,:,:],axis=1)

  # Creating netcdf file to write to
  dates = files[n][len(dirQin)+len(fileidQin):-3]
  print("Writing data for time period: "+dates)
  filename = "ERA5.MFCONV."+dates+".nc"

  print("Creating file: "+datadirout+filename)
  fileout = Dataset(datadirout+filename,"w",format="NETCDF4")

  # Global Attributes
  fileout.long_name   = "Moisture flux convergence"
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

  mfns.write_var("latitude","Latitude","","latitude",
   np.float64,"degreesN",fileout,lats)

  mfns.write_var("longitude","Longitude","","longitude",
   np.float64,"degreesE",fileout,lons)

  # Create variables in file
  mfns.write_var("CONV18","Convergence",
   "1000-850 hPa average",("time","latitude","longitude"),
   np.float64,"s**-1",fileout,conv18)

  # Create variables in file
  mfns.write_var("CONV31","Convergence",
   "300-100 hPa average",("time","latitude","longitude"),
   np.float64,"s**-1",fileout,conv31)

  # Create variables in file
  mfns.write_var("MADV18","Moisture Advection",
   "1000-850 hPa average. 1st term in MFC.",
   ("time","latitude","longitude"),
   np.float64,"kg kg**-1 s**-1",fileout,Qadv18)

  # Create variables in file
  mfns.write_var("MADV84","Moisture Advection",
   "800-400 hPa average. 1st term in MFC.",
   ("time","latitude","longitude"),
   np.float64,"kg kg**-1 s**-1",fileout,Qadv84)

  # Create variables in file
  mfns.write_var("MCONV18","Moisture Convergence",
   "1000-850 hPa average. 2nd term in MFC.",
   ("time","latitude","longitude"),
   np.float64,"kg kg**-1 s**-1",fileout,Qconv18)

  # Create variables in file
  mfns.write_var("MCONV84","Moisture Convergence",
   "800-400 hPa average. 2nd term in MFC.",
   ("time","latitude","longitude"),
   np.float64,"kg kg**-1 s**-1",fileout,Qconv84)

  # Create variables in file
  mfns.write_var("MFC18","Moisture flux convergence",
   "1000-850 hPa average",("time","latitude","longitude"),
   np.float64,"kg kg**-1 s**-1",fileout,Qfc18)

  # Create variables in file
  mfns.write_var("MFC84","Moisture flux convergence",
   "800-400 hPa average",("time","latitude","longitude"),
   np.float64,"kg kg**-1 s**-1",fileout,Qfc84)

  fileout.close()

# Loop over all desired months
#for n in range(len(files)):
#  write_calc(n)
#  exit()

Parallel(n_jobs=4)(delayed(write_calc)(n) 
  for n in range(len(files)))

