#==================================================================
# Compile IMERG and FiT information into precipitation 
# feature netcdf files. Requires driver_processFiTobs.py be 
# in same directory.
#
# James Russell 2019
#==================================================================

# Import python libraries (do not change)
import glob
import time as tm
import datetime as dt
from joblib import Parallel, delayed
import os

import sys
import numpy as np
from netCDF4 import Dataset
import scipy.ndimage as ndimage

from collections import OrderedDict as OD

import warnings
warnings.filterwarnings("ignore")

# Import namelist
from namelist_TIMPS import general as gnl
from namelist_TIMPS import process as pnl

# Load custom libraries
sys.path.insert(0,gnl.fnsdir)
import misc_functions as mfns
import earth_functions as efns

#==================================================================
# Initialize timer
#==================================================================

startnow = tm.time()

#==================================================================
# Read FiT text file
#==================================================================

# Get domain size
datasetFin = Dataset(
 f'{gnl.datadirFiTin}{gnl.fileidFiTin}00000.nc')
mxdomy = datasetFin.dimensions["y"].size-1
mxdomx = datasetFin.dimensions["x"].size-1
datasetFin.close()

# Reads directory and filenames
filenames = sorted(glob.glob(f'{gnl.datadirtrkout}/*/\
{gnl.fileidtrkout}*.nc'))
if len(filenames)==0: raise ValueError(
 "No tracking files found. Incorrect file path possible.")
lentrkstr = len(gnl.datadirtrkout)+7+len(gnl.fileidtrkout)

# Find id of first and last time
#!! Change this code for different directory structure !!#
if pnl.subsetdts:
  # Search for ID of first and last time
  for i,f in enumerate(filenames):
    if f'{f[lentrkstr:lentrkstr+12]}'==pnl.date1: 
      timeidstart = i
    if f'{f[lentrkstr:lentrkstr+12]}'==pnl.date2: 
      timeidend = i
      break
  filenums = [f[-8:-3] for f in filenames]
  timeidstart = int(filenums[timeidstart])
  timeidend = int(filenums[timeidend])

# Read text file
print("Reading FiT text file")
f = open(f'{gnl.datadirFiTout}{gnl.fileidFiTout1}\
{gnl.fileidFiToutT}','r')
FiTinfo = np.genfromtxt(f, dtype="str",
 delimiter="\t",skip_header=1)
f.close()

print("Counting objects")

# Get unique IDs and times
FiTids = FiTinfo[:,0].astype(int)
FiTtms = FiTinfo[:,2].astype(int)
obidtm = OD(sorted({i: [] for i in set(FiTids)}.items()))
for k, v in zip(FiTids, FiTtms): obidtm[k].append(v)
objs = list(obidtm.keys())

# Get first and last time of each object
obidft = OD(sorted({i: min(obidtm[i]) \
 for i in set(FiTids)}.items()))
obidlt = OD(sorted({i: max(obidtm[i]) \
 for i in set(FiTids)}.items()))

# Get number of times each object survives
obnt = list(map(len,obidtm.values()))

print("Getting extent of each object")

# Get meridional extent
FiTy = list(zip(*[map(int,i.replace('y:','').split('-')) 
 for i in FiTinfo[:,4]]))
FiTmny = FiTy[0]
FiTmxy = FiTy[1]
obidmny = OD(sorted({i: [] for i in set(FiTids)}.items()))
for k, v in zip(FiTids, FiTmny): obidmny[k].append(v)
obidmxy = OD(sorted({i: [] for i in set(FiTids)}.items()))
for k, v in zip(FiTids, FiTmxy): obidmxy[k].append(v)

# Get zonal extent
FiTx = list(zip(*[map(int,i.replace('x:','').split('-')) 
 for i in FiTinfo[:,3]]))
FiTmnx = FiTx[0]
FiTmxx = FiTx[1]
obidmnx = OD(sorted({i: [] for i in set(FiTids)}.items()))
for k, v in zip(FiTids, FiTmnx): obidmnx[k].append(v)
obidmxx = OD(sorted({i: [] for i in set(FiTids)}.items()))
for k, v in zip(FiTids, FiTmxx): obidmxx[k].append(v)

# Find only obids of objects that reach a size of nsz 
#  pixels in their lifetime
if pnl.subsetsz:
  print("Subsetting by size")

  # Sort sizes
  FiTszs = [int(i.split(":")[1]) for i in FiTinfo[:,8]]
  obidsz = OD(sorted({i: [] for i in set(FiTids)}.items()))
  for k, v in zip(FiTids, FiTszs): obidsz[k].append(v)

  obmsz = list(map(max,obidsz.values()))
  obidmsz = OD(sorted({i: [] for i in set(FiTids)}.items()))
  for k, v in zip(objs, obmsz): obidmsz[k] = v
 
  # Find indices where larger than size
  indgtmsz = np.where(np.array(obmsz)>=pnl.nsz)[0]
  #keyltmsz = [k for k, v in obidmsz.items() if v<pnl.nsz]

  # Keep only those objects that are larger
  objs = [objs[i] for i in indgtmsz]
  obidft = {i: obidft[i] for i in objs}
  obidlt = {i: obidlt[i] for i in objs}
  obnt = [obnt[i] for i in indgtmsz]
  obidmxx = {i: obidmxx[i] for i in objs}
  obidmnx = {i: obidmnx[i] for i in objs}
  obidmxy = {i: obidmxy[i] for i in objs}
  obidmny = {i: obidmny[i] for i in objs}

# Find only obids of objects that exist for at least nt
#  times.
if pnl.subsettm:
  print("Subsetting by time")

  # Find indices where larger than size
  indgtmtm = np.where(np.array(obnt)>=pnl.nt)[0]

  # Keep only those objects that last longer
  objs = [objs[i] for i in indgtmtm]
  obidft = {i: obidft[i] for i in objs}
  obidlt = {i: obidlt[i] for i in objs}
  obidmxx = {i: obidmxx[i] for i in objs}
  obidmnx = {i: obidmnx[i] for i in objs}
  obidmxy = {i: obidmxy[i] for i in objs}
  obidmny = {i: obidmny[i] for i in objs}

# Subset by range of object ids
if pnl.subsetobs:
  print("Subsetting by objects")
  objs = list(np.intersect1d(
              np.array(np.arange(pnl.ob1,pnl.ob2)),
              np.array(objs)))
  obidft = {i: obidft[i] for i in objs}
  obidlt = {i: obidlt[i] for i in objs}
  obidmxx = {i: obidmxx[i] for i in objs}
  obidmnx = {i: obidmnx[i] for i in objs}
  obidmxy = {i: obidmxy[i] for i in objs}
  obidmny = {i: obidmny[i] for i in objs}

# Find only obids of objects that are not on domain boundary

# Find objects and indices
oldobjs = objs
objsnotonmny = [k for k,v in obidmny.items() if min(v)>0]
objsnotonmxy = [k for k,v in obidmxy.items() if max(v)<mxdomy]
objs = np.intersect1d(objsnotonmny,objsnotonmxy)
if not gnl.domperiodicx:
  objsnotonmnx = [k for k,v in obidmnx.items() if min(v)>0]
  objsnotonmxx = [k for k,v in obidmxx.items() if max(v)<mxdomx]
  objsx = np.intersect1d(objsnotonmnx,objsnotonmnx)
  objs = np.intersect1d(objs,objsx)

# Keep only those objects that are not on domain boundary
indnotondb = np.intersect1d(oldobjs,objs,return_indices=True)[1]
obnt = [obnt[i] for i in indnotondb]
obidft = {i: obidft[i] for i in objs}
obidlt = {i: obidlt[i] for i in objs}
obidmxx = {i: obidmxx[i] for i in objs}
obidmnx = {i: obidmnx[i] for i in objs}
obidmxy = {i: obidmxy[i] for i in objs}
obidmny = {i: obidmny[i] for i in objs}

# Find only obids of objects that start between date1 and date2

# Find objects and indices
if pnl.subsetdts:
  oldobjs = objs
  objs = [i for i,v in obidft.items() if timeidstart<=v<=timeidend]

  # Keep only those objects that are not on domain boundary
  indnotondb = np.intersect1d(oldobjs,objs,return_indices=True)[1]
  obnt = [obnt[i] for i in indnotondb]
  obidft = {i: obidft[i] for i in objs}
  obidlt = {i: obidlt[i] for i in objs}
  obidmxx = {i: obidmxx[i] for i in objs}
  obidmnx = {i: obidmnx[i] for i in objs}
  obidmxy = {i: obidmxy[i] for i in objs}
  obidmny = {i: obidmny[i] for i in objs}

print("Remove all objects with first time during spinup")
# Find indices of objects starting before spinup
indspin = np.where(np.array(list(obidft.values()))>=int(filenames[0][-8:-3]))[0]

# Keep only those objects that are larger
objs = [objs[i] for i in indspin]
obidft = {i: obidft[i] for i in objs}
obidlt = {i: obidlt[i] for i in objs}
obnt = [obnt[i] for i in indspin]
obidmxx = {i: obidmxx[i] for i in objs}
obidmnx = {i: obidmnx[i] for i in objs}
obidmxy = {i: obidmxy[i] for i in objs}
obidmny = {i: obidmny[i] for i in objs}

print(f'Removed {str(len(obidtm.keys())-len(objs))}\
 of {str(len(obidtm.keys()))} objects. Processing\
 {len(objs)} objects.')

del(FiTinfo)

# Obtain universal information
datasetF0 = Dataset(filenames[0])
lat = datasetF0.variables["lat"][:]
lon = datasetF0.variables["lon"][:]
datasetF0.close()

# Set reference time
# Set reference time
sincedate = dt.datetime.strptime(pnl.reftime,
 "%Y-%m-%d %H:%M:%S")

#==================================================================
# Define function
#==================================================================

def driver_processFiTobs(o):
  "o corresponds to the object in objs"

  warnings.filterwarnings("ignore")

#==================================================================
# Initialize function variables
#==================================================================

  times = np.arange(obidft[objs[o]],obidlt[objs[o]]+1,1)
  print(f'Working on object {str(objs[o])} with\
 {str(len(times))} times')

  # Times in hours since reftime
  time1          = [np.nan]* len(times) 
  # Times (YYYYMMDDhhmm)
  datetime       = [np.nan] * len(times)
  # Central latitude of object
  centrallat     = [np.nan] * len(times) 
  # Central longitude of object
  centrallon     = [np.nan] * len(times)
  # Weighted central latitude of object
  wgtcentlat     = [np.nan] * len(times)
  # Weighted central longitude of object
  wgtcentlon     = [np.nan] * len(times)

  # Dictionaries
  lats = {}; lons = {}; instrain = {}#; QI = {}

  # Loop over all times relevant to object
  for t in range(0,len(times)):

#==================================================================
# Read in data for current FiT objects

    # Convert latitudes
    latnow = lat[obidmny[objs[o]][t]:obidmxy[objs[o]][t]+1]
    lonnow = lon[obidmnx[objs[o]][t]:obidmxx[objs[o]][t]+1]

    # Open file and read map of object ids
    trkfilename = glob.glob(f'{gnl.datadirtrkout}*/\
{gnl.fileidtrkout}*_{str(times[t]).zfill(5)}.nc')[0]
    timenow    = trkfilename[lentrkstr:lentrkstr+12]
    datasettrk = Dataset(trkfilename)
    objects    = np.squeeze(datasettrk.variables["SystemID"][:,
     obidmny[objs[o]][t]:obidmxy[objs[o]][t]+1,
     obidmnx[objs[o]][t]:obidmxx[objs[o]][t]+1])

#==================================================================
# Get indices for current object

    # If 1D, add singleton dimension
    if obidmny[objs[o]][t]==obidmxy[objs[o]][t]:
      objects = np.expand_dims(objects,0)
    if obidmnx[objs[o]][t]==obidmxx[objs[o]][t]:
      objects = np.expand_dims(objects,1)

    # Find all location indices of object
    locindob = np.where(np.where(
     objects==int(objs[o]),1,0)==1)

#==================================================================
# Read IMERG data and process for specific object

    # Read variables
    rain = datasettrk.variables["Precipitation"][:,
     obidmny[objs[o]][t]:obidmxy[objs[o]][t]+1,
     obidmnx[objs[o]][t]:obidmxx[objs[o]][t]+1]
    datasettrk.close()

    # Change dimension order of variable
    rain = np.squeeze(rain)

    # If 1D, add singleton dimension
    if obidmny[objs[o]][t]==obidmxy[objs[o]][t]:
      rain = np.expand_dims(rain,0)
    if obidmnx[objs[o]][t]==obidmxx[objs[o]][t]:
      rain = np.expand_dims(rain,1)

#==================================================================
# Get rain, lat, and lon of current objects
   
    rain1 = [rain[y,x] for y,x in 
     list(zip(locindob[0],locindob[1]))]
    lat1 = [latnow[l] for l in locindob[0]]
    lon1 = [lonnow[l] for l in locindob[1]]
   
#==================================================================
# Calculate parameters and assign data at first time
   
    # Assign some simple values
    datetime[t] = int(timenow)

    # Calculate center of mass in lat
    centrallat[t]  = np.mean(lat1)
    
    # Calculated differently for periodic boundaries in lon
    centrallon[t] = efns.periodic_cmass(lon1)

    # Calculate weighted (by rainfall) centroid
    sr = sum(rain1)
    if sr!=0:
      rw = [mfns.divzero(r,sr) for r in rain1]
      wgtcentlat[t] = np.average(lat1,weights=rw)
      wgtcentlon[t] = efns.periodic_cmass_weighted(lon1,rw)

    # Calculate time in hours since reftime
    currentdate = dt.datetime.strptime(timenow,"%Y%m%d%H%M")
    delta = (currentdate - sincedate)
    time1[t] = delta.total_seconds()/3600.

    # For each new time, add to dictionary as a new key
    lats[str(datetime[t])] = lat1
    lons[str(datetime[t])] = lon1
    instrain[str(datetime[t])] = rain1

  # If no precipitation value passes the threshold, ignore
  if pnl.subsetmrn:
    passrthold = False
    for val in instrain.values():
      if np.amax(val)>gnl.convrainthold:
        passrthold = True
        continue
    if not passrthold:
      print("System does not reach minimum rain threshold")
      return

  # If no time passes the minimum size, ignore
  if pnl.subsetsz:
    passsthold = False
    for val in instrain.values():
      val = np.array(val)
      if len(val[val>0])>pnl.nsz:
        passsthold = True
        continue
    if not passsthold:
      print("System does not reach minimum size")
      return

  # For PF, if all rain rates at all times are actually zero
  #  move onto the next PF
  for val in instrain.values():
    if all(v == 0 for v in val):
      a = True
    else:
      a = False
      break
  if a:
    return
  
  # For all other PFs, if all rain rates at any given time 
  #  are zero, add that index to a list of indices to delete
  delind = []
  delkeys = []
  for key,val in instrain.items():
    if all(v == 0 for v in val):
      delkeys.append(key)
  setdk = set(delkeys)
  delind = [i for i, e in enumerate(instrain.keys()) 
   if e in setdk]

  # Delete the key, value pairs in the dictionaries
  for i in delind:
    lats.pop(str(datetime[i]))
    lons.pop(str(datetime[i]))
    instrain.pop(str(datetime[i]))

  # If not enough times
  if pnl.subsettm:
    if len(instrain)<pnl.nt:
      print("System does not reach minimum # times")
      return

  # Delete the indices in the lists
  delind.reverse()
  datetime = [i for j, i in enumerate(datetime) 
   if j not in delind]
  time1 = [i for j, i in enumerate(time1) 
   if j not in delind]
  centrallat = [i for j, i in enumerate(centrallat) 
   if j not in delind]
  centrallon = [i for j, i in enumerate(centrallon) 
   if j not in delind]
  wgtcentlat = [i for j, i in enumerate(wgtcentlat) 
   if j not in delind]
  wgtcentlon = [i for j, i in enumerate(wgtcentlon) 
   if j not in delind]

  # Convert lons to line up with rest
  centrallon = [l-360 if l>180 else l for l in centrallon]
  wgtcentlon = [l-360 if l>180 else l for l in wgtcentlon]

#==================================================================
# Write NetCDF file for object

  # Creating netcdf file
  datadir = f'{gnl.datadirTIMPS}{str(datetime[0])[4:6].zfill(2)}/'
  try: os.mkdir(datadir)
  except: print("Directory already exists")
  filename = f'{gnl.fileidTIMPS}{str(objs[o]).zfill(7)}_\
{str(datetime[0]).zfill(4)}_\
{str(int(round(centrallat[0]))).zfill(2)}_\
{str(int(round(centrallon[0]))).zfill(2)}.nc'

  print(f'Creating file: {datadir}{filename}')
  fileout = Dataset(f'{datadir}{filename}',
   'w',format='NETCDF4')

  # Global Attributes
  fileout.long_name   = 'Precipitation feature information'
  fileout.description = 'Information for a precipitation\
 feature tracked by FiT algorithm with specifications as\
 in the attributes below'
  fileout.history     = f'Created {str(tm.ctime(tm.time()))}'
  fileout.SystemID    = str(objs[o])

  # Create dimensions in file
  time = fileout.createDimension('time', len(time1))

  # Create variables in file
  mfns.write_var("time","Time","","time",np.float64,
   f'hours since {pnl.reftime}',fileout,time1)

  mfns.write_var("datetime","Date and time","","time",
   np.int64,"YYYYMMDDhhmm",fileout,datetime)

  mfns.write_var("centrallat","Central latitude",
   "Latitude of PF centroid","time",np.float64,
   "degreesNorth",fileout,centrallat)

  mfns.write_var("centrallon","Central longitude",
   "Longitude of PF centroid","time",np.float64,
   "degreesEast",fileout,centrallon)

  description = \
   "Latitude of PF centroid weighted by rainfall"
  mfns.write_var("centlatwgt","Weighted Central Latitude",
   description,"time",np.float64,"degreesNorth",
   fileout,wgtcentlat)

  description = \
   "Longitude of PF centroid weighted by rainfall"
  mfns.write_var("centlonwgt","Weighted Central Longitude",
   description,"time",np.float64,"degreesEast",
   fileout,wgtcentlon)

  format1 = "Data is in attribute and value pairs of the\
 subgroup data. Attributes correspond to the date and time\
 in YYYYMMDDhhmm format. Values of those attributes are\
 lists of the data at that time."

  mfns.write_group("lats","Latitudes",
   "Latitudes of IMERG grid cell centers for PF",
   "DegreesNorth",format1,fileout,
   lats,f'{gnl.datadirTIMPS}{filename}')

  mfns.write_group("lons","Longitudes",
   "Longitudes of IMERG grid cell centers for PF",
   "DegreesEast",format1,fileout,
   lons,f'{gnl.datadirTIMPS}{filename}')

  description = "Instantaneous rain rates at IMERG grid\
cells corresponding to all latitude and longitude values in PF."
  mfns.write_group("instrain","Instantaneous rain rate",
                 description,"mm/hr",format1,fileout,
                 instrain,f'{gnl.datadirTIMPS}{filename}')
  # Close file
  fileout.close()

#==================================================================
# End Function
#================================================================== 

#==================================================================
# Parallel loop over object
#================================================================== 

# Begins loop
if pnl.serialorparallel==1:
  print(f"Begin serial loop over  {len(objs)} objects")
  for o in range(len(objs)):
    driver_processFiTobs(o)

# Parallel loop over PFs
if pnl.serialorparallel==2:
  print(f"Begin parallel loop over {len(objs)} objects")
  Parallel(n_jobs=int(pnl.njobs))(delayed(
   driver_processFiTobs)(o) for o in range(len(objs)))

#==================================================================
# Final clean up tasks
#==================================================================

# Remove namelist and filename files
print("Cleaning up")
os.system("rm -rf __pycache__")

# End timer
endnow = tm.time()
print(f'Program took {str(endnow - startnow)}s')

#==================================================================
# End code
#==================================================================
