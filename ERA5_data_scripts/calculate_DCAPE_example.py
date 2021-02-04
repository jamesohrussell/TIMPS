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
import sharppy.sharptab.params  as params

# Namelist
datadirout = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/CAPE/"

dirTin     = "/uufs/chpc.utah.edu/common/home/u0816744/ERA5/temperature/allplevs/"
fileidTin  = "ERA5.TEMP."
dirRin     = "/uufs/chpc.utah.edu/common/home/u0816744/ERA5/moisture/RLHU/allplevs/"
fileidRin  = "ERA5.RLHU."
dirGin     = "/uufs/chpc.utah.edu/common/home/u0816744/ERA5/geopotential/allplevs/"
fileidGin  = "ERA5.GEOP."

latS = -35
latN = 35
levh = 100
levl = 1000

# Find all temperature files
print("Finding all files to process")
Tfiles = sorted(glob.glob(dirTin+fileidTin+"*"))

def write_calc(n):

  print("Working on file "+Tfiles[n])

  Tf = Dataset(Tfiles[n])

  # Get pressure and coordinates
  print("Getting pressure and coordinates")
  time = Tf.variables["time"][:]
  levs = Tf.variables["level"][:]
  lats = Tf.variables["latitude"][:]
  lons = Tf.variables["longitude"][:]

  latiN = mfns.k_closest(lats,latN,1)[0]
  latiS = mfns.k_closest(lats,latS,1)[0]
  levih = mfns.k_closest(levs,levh,1)[0]
  levil = mfns.k_closest(levs,levl,1)[0]

  P = levs[levih:levil+1]
  lats = lats[latiN:latiS+1]

  # Get temperature
  print("Getting temperature")
  T = Tf.variables["T"][:,levih:levil+1,latiN:latiS+1,:]-273.15

  # Get dew point
  print("Getting dewpoint")
  Rf = Dataset(dirRin+fileidRin+\
   Tfiles[n][len(dirTin)+len(fileidTin):-3]+".nc")
  R = Rf.variables["R"][:,levih:levil+1,latiN:latiS+1,:]
  alpha = ((17.27*T)/(237.7+T))+np.log(R/100.)
  D = (237.7*alpha)/(17.27-alpha)

  # Get height
  print("Getting height")
  Gf = Dataset(dirGin+fileidGin+\
   Tfiles[n][len(dirTin)+len(fileidTin):-3]+".nc")
  Z = Gf.variables["Z"][:,levih:levil+1,latiN:latiS+1,:]/9.80665

  # Make bogus wind speed and direction profiles
  WS = np.zeros(P.shape)
  WD = np.zeros(P.shape)

  # Calculate DCAPE
  print("Calculating DCAPE")
  P = P[::-1]
  T = T[:,::-1,:,:]
  D = D[:,::-1,:,:]
  Z = Z[:,::-1,:,:]

  print(P)
  print(T[0,:,0,0])
  print(D[0,:,0,0])
  print(Z[0,:,0,0])

  Tshp = T.shape
  DCAPE = np.array([[[params.dcape(profile.create_profile(profile='default',
   pres=P,hght=Z[t,:,j,i],tmpc=T[t,:,j,i],dwpc=D[t,:,j,i],wspd=WS,wdir=WD))[0]
   for t in range(Tshp[0])] for j in range(Tshp[2])] for i in range(Tshp[3])])

  print(DCAPE.shape)
  exit()

  #params.dcape(profile.create_profile(profile='default',
  # pres=P,hght=Z,tmpc=T,dwpc=D,wspd=WS,wdir=WD))[0]


Parallel(n_jobs=1)(delayed(write_calc)(n) 
  for n in range(len(Tfiles)))

#alpha = ((17.27*T-273.15)/(237.7+T-273.15))+np.log(R/100.)
#D = (237.7*alpha)/(17.27-alpha)
#
#
#P  = [1000,900,800,700,600,500,400,300,200,100]
#Z  = [0,1000,2000,3000,4500,6000,8000,10000,13000,17000]
#T  = [27,22,17,13,5,-3,-13,-28,-50,-80]
#D  = [23,10,5, -5,-15,-20,-35,-50,-75,-110]
#WS = [0,0,0,0,0,0,0,0,0,0]
#WD = [0,0,0,0,0,0,0,0,0,0]
#
#dcape = params.dcape(profile.create_profile(profile='default',
# pres=P,hght=Z,tmpc=T,dwpc=D,wspd=WS,wdir=WD))[0]
#
#print(dcape)
#exit()
