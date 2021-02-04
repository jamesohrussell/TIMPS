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
import h5py

from collections import OrderedDict as OD

import warnings
warnings.filterwarnings("ignore")

# Import namelist
from namelist_TIMPS import general as gnl
from namelist_TIMPS import process as pnl

# Load custom libraries
sys.path.insert(0,gnl.fnsdir)
import misc_functions as mfns
import shape_functions as sfns

#==================================================================
# Initialize timer
#==================================================================

startnow = tm.time()

#==================================================================
# Get universal information
#==================================================================

# Get all original files used for tracking
print("Generating list of files")
datasetFin = Dataset(
 f'{gnl.datadirFiTin}{gnl.fileidFiTin}00000.nc')
datestart = datasetFin.datestart
dateend = datasetFin.dateend
mxdomy = datasetFin.dimensions["y"].size-1
mxdomx = datasetFin.dimensions["x"].size-1

# Generate a list of filenames with dates to search for
start = dt.datetime.strptime(datestart, "%Y%m%d")
end = dt.datetime.strptime(dateend, "%Y%m%d")
datearr = (start + dt.timedelta(days=x) \
          for x in range(0, (end-start).days))
filen = []
#!! Change this line for different directory structure !!#
for dateobj in datearr:
  filen.append(f'{gnl.datadirin}\
{dateobj.strftime("%Y")}/{dateobj.strftime("%m")}/\
{gnl.fileidin+dateobj.strftime("%Y%m%d")}*')

# Reads directory and filenames
n = 0
for f in filen:
  if n==0:
    filenames = glob.glob(f)
  else:
    filenames = filenames+glob.glob(f)
  n=n+1
lenIdir = len(gnl.datadirin)

# Find id of first and last time
#!! Change this code for different directory structure !!#
lenstr = len(f'{gnl.datadirin}{dateobj.strftime("%Y")}/\
{dateobj.strftime("%m")}/{gnl.fileidin}')
if pnl.subsetdts:
  # Search for ID of first and last time
  for i,f in enumerate(filenames):
    if f'{f[lenstr:lenstr+8]}\
{f[lenstr+10:lenstr+14]}'==pnl.date1:
      timeidstart = i
    if f'{f[lenstr:lenstr+8]}\
{f[lenstr+10:lenstr+14]}'==pnl.date2:
      timeidend = i
      break

#==================================================================
# Read FiT text file
#==================================================================

# Read text file
print("Reading FiT text file")
f = open(f'{gnl.datadirFiTout}{gnl.fileidFiTout1}\
{gnl.fileidFiToutT}','r')
FiTinfo = np.genfromtxt(f, dtype="str",delimiter="\t",
                           skip_header=1)
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
  objs = [i for i,v in obidft.items() 
   if timeidstart<=v<=timeidend]

  # Keep only those objects that are not on domain boundary
  indnotondb = np.intersect1d(oldobjs,objs,
   return_indices=True)[1]
  obnt = [obnt[i] for i in indnotondb]
  obidft = {i: obidft[i] for i in objs}
  obidlt = {i: obidlt[i] for i in objs}
  obidmxx = {i: obidmxx[i] for i in objs}
  obidmnx = {i: obidmnx[i] for i in objs}
  obidmxy = {i: obidmxy[i] for i in objs}
  obidmny = {i: obidmny[i] for i in objs}

print(f'Removed {str(len(obidtm.keys())-len(objs))}\
 of {str(len(obidtm.keys()))} objects. Processing\
 {len(objs)} objects.')

del(FiTinfo)

# Get universal information (lat,lon,subset info)
finfo = Dataset(f'{gnl.datadirFiTin}\
{gnl.fileidFiTin}00000.nc')

# Obtain universal information
region = datasetFin.region
tholds = datasetFin.tholds
source = datasetFin.source
if region=='global':
  ssreg    = False
else:
  ssreg    = True
  latN     = datasetFin.latN
  latS     = datasetFin.latS
  lonW     = datasetFin.lonW
  lonE     = datasetFin.lonE
if hasattr(datasetFin, 'smthstddev'):
  sstd     = datasetFin.smthstddev
else:
  sstd     = None
if hasattr(datasetFin, 'smthwidth'):
  swid     = datasetFin.smthwidth
else:
  swid     = None

datasetFin.close()

if region!="global":

  # Read first file
  datasetIM = h5py.File(filenames[0],'r')

  # Read in coordinate variables
  lat = datasetIM["Grid/lat"][:]
  lon = datasetIM["Grid/lon"][:]
  datasetIM.close()
        
  # Find indices of closest coordinates
  idyS = (np.abs(lat - latS)).argmin()
  idyN = (np.abs(lat - latN)).argmin()
  idxW = (np.abs(lon - lonW)).argmin()
  idxE = (np.abs(lon - lonE)).argmin()

  # Subset main dataset by indices
  lat  = lat[idyS:idyN+1]
  lon  = lon[idxW:idxE+1]

finfo.close()

# Set reference time
sincedate = dt.datetime.strptime(pnl.reftime,
 "%Y-%m-%d %H:%M:%S")

#==================================================================
# Define function
#==================================================================

def driver_processFiTobs(o):
  "o corresponds to the object in objs"

#==================================================================
# Initialize function variables
#==================================================================

  times = np.arange(obidft[objs[o]],obidlt[objs[o]]+1,1)
  print(f'Working on object {str(objs[o])} with\
 {str(len(times))} times')

  # Times in hours since reftime
  time1          = [np.nan]*len(times)
  # Times (YYYYMMDDhhmm)
  datetime       = [np.nan]*len(times)
  # Central latitude of object
  centrallat     = [np.nan]*len(times)
  # Central longitude of object
  centrallon     = [np.nan]*len(times)
  # Weighted central latitude of object
  wgtcentlat     = [np.nan]*len(times)
  # Weighted central longitude of object
  wgtcentlon     = [np.nan]*len(times)

  # Dictionaries
  lats = {}; lons = {}; instrain = {}#; QI = {}

  # Loop over all times relevant to object
  for t in range(0,len(times)):

#==================================================================
# Read in data for current FiT objects

    # Open FiT input file
    Finfilename = f'{gnl.datadirFiTin}{gnl.fileidFiTin}\
{str(times[t]).zfill(5)}.nc'
    datasetFin  = Dataset(Finfilename)
    timenow     = datasetFin.time
    datasetFin.close()

    # Convert latitudes
    latnow = lat[obidmny[objs[o]][t]:obidmxy[objs[o]][t]+1]
    lonnow = lon[obidmnx[objs[o]][t]:obidmxx[objs[o]][t]+1]

    # Open file and read map of object ids
    FiTfilename = f'{gnl.datadirFiTout}{gnl.fileidFiTout1}\
{str(times[t]).zfill(5)}{gnl.fileidFiTout2}'
    datasetFi   = Dataset(FiTfilename)
    objects = np.squeeze(datasetFi.variables["value"][:,
     obidmny[objs[o]][t]:obidmxy[objs[o]][t]+1,
     obidmnx[objs[o]][t]:obidmxx[objs[o]][t]+1])
    datasetFi.close()

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
  
    # Open file
    datasetIM = h5py.File(filenames[times[t]],'r')
    # Read variables
    rain = datasetIM["/Grid/precipitationCal"][:,
     idxW+obidmnx[objs[o]][t]:idxW+obidmxx[objs[o]][t]+1,
     idyS+obidmny[objs[o]][t]:idyS+obidmxy[objs[o]][t]+1]
    datasetIM.close()

    # Change dimension order of variable
    rain = np.squeeze(np.transpose(rain,(0, 2, 1)))

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
    xiibar = np.cos(np.radians(lon1))
    zetaibar = np.sin(np.radians(lon1))

    # Unweighted central lon
    centrallon[t] = np.degrees(np.pi+
     np.arctan2(-np.mean(zetaibar),-np.mean(xiibar)))

    # Calculate weighted (by rainfall) centroid
    sr = sum(rain1)
    if sr!=0:
      rw = [mfns.divzero(r,sr) for r in rain1]
      wgtcentlat[t] = np.average(lat1,weights=rw)
      wgtcentlon[t] = np.degrees(np.pi+
         np.arctan2(-np.average(zetaibar,weights=rw),
                    -np.average(xiibar,weights=rw)))

    # Calculate time in hours since reftime
    currentdate = dt.datetime.strptime(timenow,"%Y%m%d%H%M")
    delta = (currentdate - sincedate)
    time1[t] = delta.total_seconds()/3600.

    # For each new time, add to dictionary as a new key
    lats[str(datetime[t])] = lat1
    lons[str(datetime[t])] = lon1
    instrain[str(datetime[t])] = rain1

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

  # Creating netcdf file to write to
  filename = f'{gnl.fileidTIPS}{str(objs[o]).zfill(7)}_\
{str(datetime[0])}_\
{str(int(round(centrallat[0]))).zfill(2)}_\
{str(int(round(centrallon[0]))).zfill(2)}.nc'

  print(f'Creating file: {gnl.datadirTIPS}{filename}')
  fileout = Dataset(f'{gnl.datadirTIPS}{filename}',
   'w',format='NETCDF4')

  # Global Attributes
  fileout.long_name   = 'Precipitation feature information'
  fileout.description = 'Information for a precipitation\
 feature tracked by FiT algorithm'
  fileout.history     = f'Created {str(tm.ctime(tm.time()))}'
  fileout.systemID    = str(objs[o])

  # Create dimensions in file
  time = fileout.createDimension('time', len(time1))

  # Create variables in file
  mfns.write_var("time","Time","","time",np.float64,
   f'hours since {pnl.reftime}',fileout,time1,
   f'{gnl.datadirTIPS}{filename}')

  mfns.write_var("datetime","Date and time","","time",
   np.int64,"YYYYMMDDhhmm",fileout,datetime,
   f'{gnl.datadirTIPS}{filename}')

  mfns.write_var("centrallat","Central latitude",
   "Latitude of PF centroid","time",np.float64,
   "degreesNorth",fileout,centrallat,
   f'{gnl.datadirTIPS}{filename}')

  mfns.write_var("centrallon","Central longitude",
   "Longitude of PF centroid","time",np.float64,
   "degreesEast",fileout,centrallon,
   f'{gnl.datadirTIPS}{filename}')

  description = \
   "Latitude of PF centroid weighted by rainfall"
  mfns.write_var("centlatwgt","Weighted Central Latitude",
   description,"time",np.float64,"degreesNorth",fileout,
   wgtcentlat,f'{gnl.datadirTIPS}{filename}')

  description = \
   "Longitude of PF centroid weighted by rainfall"
  mfns.write_var("centlonwgt","Weighted Central Longitude",
   description,"time",np.float64,"degreesEast",fileout,
   wgtcentlon,f'{gnl.datadirTIPS}{filename}')

  format1 = "Data is in attribute and value pairs of the\
 subgroup data. Attributes correspond to the date and time\
 in YYYYMMDDhhmm format. Values of those attributes are\
 lists of the data at that time."

  mfns.write_group("lats","Latitudes",
   "Latitudes of IMERG grid cell centers for PF",
   "DegreesNorth",format1,fileout,
   lats,f'{gnl.datadirTIPS}{filename}')

  mfns.write_group("lons","Longitudes",
   "Longitudes of IMERG grid cell centers for PF",
   "DegreesEast",format1,fileout,
   lons,f'{gnl.datadirTIPS}{filename}')

  description = "Instantaneous rain rates at IMERG grid\
cells corresponding to all latitude and longitude values in PF."
  mfns.write_group("instrain","Instantaneous rain rate",
                 description,"mm/hr",format1,fileout,
                 instrain,f'{gnl.datadirTIPS}{filename}')
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
  print("Begin serial loop over objects")
  for o in range(len(objs)):
    driver_processFiTobs(o)

# Parallel loop over PFs
if pnl.serialorparallel==2:
  print("Begin parallel loop over objects")
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
