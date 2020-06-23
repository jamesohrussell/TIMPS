#============================================================
# Compile IMERG and FiT information into precipitation 
# feature netcdf files. Requires driver_processFiTobs.py be 
# in same directory.
#
# James Russell 2019
#============================================================

# Import python libraries (do not change)
import numpy as np
from netCDF4 import Dataset
import glob
from collections import Counter
import datetime as dt
import driver_processFiTobs as dp
from joblib import Parallel, delayed
import time as tm
import os

# Import namelist
from namelist_TIPS import general as gnl
from namelist_TIPS import process as pnl

#============================================================
# Initialize timer
#============================================================

startnow = tm.time()

#============================================================
# Get universal information
#============================================================

# Get all original files used for tracking
print("Generating list of files")
datasetFin = Dataset(pnl.datadirFin+pnl.fileidFin+"00000.nc")
datestart = datasetFin.datestart
dateend = datasetFin.dateend

# Generate a list of filenames with dates to search for
start = dt.datetime.strptime(datestart, "%Y%m%d")
end = dt.datetime.strptime(dateend, "%Y%m%d")
datearr = (start + dt.timedelta(days=x) \
          for x in range(0, (end-start).days))
filen = []
for dateobj in datearr:
  filen.append(gnl.datadirin+gnl.fileidin+ \
  dateobj.strftime("%Y%m%d")+"*")

# Reads directory and filenames
n = 0
for f in filen:
  if n==0:
    filenames = glob.glob(f)
  else:
    filenames = filenames+glob.glob(f)
  n=n+1
lenIdir    = len(gnl.datadirin)

#============================================================
# Read FiT text file
#============================================================

# Read text file
f       = open(pnl.datadirFi+pnl.fileidFi1+pnl.fileidtxt,'r') 
FiTinfo = np.genfromtxt(f, dtype="str",delimiter="\t",
                           skip_header=1)

# Get list of unique object ids
FiTids = FiTinfo[:,0].astype(int)
FiTszs = [int(i.split(":")[1]) for i in FiTinfo[:,8]]
obids = Counter(FiTids)
objs  = list(obids.keys())

# Find only obids of objects that reach a size of nsz 
#  pixels in their lifetime and exist for at least nt
#  times.
if pnl.subsetsztm:
  obinds = [list(np.where(FiTids==key)[0]) 
           for key in objs]
  obnt  = list(obids.values())
  obmsz = [max([FiTszs[i] for i in inds]) 
                          for inds in obinds]
  obinds1 = np.intersect1d(
            np.where(np.array(obmsz)>=pnl.nsz)[0],
            np.where(np.array(obnt)>=pnl.nt)[0])
  objs = [objs[i] for i in obinds1]

# Subset by range of object ids
if pnl.subsetobs:
  objs = list(np.intersect1d(
              np.array(np.arange(pnl.ob1,pnl.ob2)),
              np.array(objs)))

#============================================================
# Write information
#============================================================ 

# Write namelist dictionary to netcdf file for reading 
#  during parallel loop
print("Writing object id to netcdf file")
pifileout = Dataset("passinfo_po.nc","w",format="NETCDF4")
pifileout.objs = objs

# Write file names to a text file for reading during 
#  parallel loop
print("Writing filenames to text file")
ffn = open("filenames_po.txt","w")
for i in range(len(filenames)):
  if i<len(filenames)-1:
    ffn.write(filenames[i]+",")
  else:
    ffn.write(filenames[i])
ffn.close()

pifileout.close()

#============================================================
# Parallel loop over object
#============================================================ 

if pnl.serialorparallel==1:
  print("Begin serial loop over objects")
  for o in range(len(objs)):
    dp.driver_processFiTobs(o)

# Parrallel loop over PFs
if pnl.serialorparallel==2:
  print("Begin parrallel loop over objects")
  Parallel(n_jobs=pnl.njobs)(delayed(dp.driver_processFiTobs)(o) \
    for o in range(len(objs)))

#============================================================
# Final clean up tasks
#============================================================

# Remove namelist and filename files
print("Cleaning up")
os.system("rm filenames_po.txt")
os.system("rm passinfo_po.nc")
os.system("rm -rf __pycache__")

# End timer
endnow = tm.time()
print("Program took "+str(endnow - startnow)+"s")

#============================================================
# End code
#============================================================
