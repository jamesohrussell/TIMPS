#============================================================
# Compile IMERG and FiT information into precipitation 
# feature netcdf files.
#
# Must provide in namelist:
#  * Directory and filename identifiers for IMERG files 
#      and FiT output netcdf files (the latter has no 
#      effect if only plotting IMERG)
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

#============================================================
# Namelist
#============================================================

# Directory for custom functions
fnsdir = "/uufs/chpc.utah.edu/common/home/u0816744/general_functions/"

# Directory and filename for input IMERG data
# Should obtain from create_FiT_input_files.py
datadirIM  = "/uufs/chpc.utah.edu/common/home/varble-group2/IMERG/"
fileidIM   = "3B-HHR.MS.MRG.3IMERG."
datahdf5   = True
datanc4    = False

# Directory and filename for FiT input files
# Should obtain from create_FiT_input_files.py
datadirFin = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_CPEX-AW/FiT_input_test/"
fileidFin  = "IMERG_FiT_tholds_CapVer_"

# Directory and filename for FiT output data files
datadirFi  = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_CPEX-AW/FiT_output_test/"
fileidFi1  = "IMERG_tracked_"
fileidFi2  = "_4Dobjects.nc"
fileidtxt  = "4Dobject_tree.txt"

# Directory and filename for output PF nc files
datadirout = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_CPEX-AW/TIPS_test/"
fileidout1 = "TIPS_"

# Only process a set of objects
subsetobs = True
ob1       = 100000 # First obid
ob2       = 100100 # Last obid
subsetsztm= True
nsz       = 10     # Minimum no of pixels
nt        = 6      # Minimum number of times

# Specify merging distance from tracking (number of pixels)
mergedist = 30

# Specify reference time in format YYYY-MM-DD hh:mm:ss
reftime = "1900-01-01 00:00:00"

# Number of processes for parallelization
serialorparallel = 2
njobs = 16

#============================================================
# Initialize timer
#============================================================

startnow = tm.time()

#============================================================
# Write namelist to a dictionary
#============================================================

# Put namelist information in dictionaries
namelist = {}
namelist["fnsdir"] = str(fnsdir)
namelist["datahdf5"] = str(datahdf5)
namelist["datanc4"] = str(datanc4)
namelist["datadirFin"] = str(datadirFin)
namelist["fileidFin"] = str(fileidFin)
namelist["datadirFi"] = str(datadirFi)
namelist["fileidFi1"] = str(fileidFi1)
namelist["fileidFi2"] = str(fileidFi2)
namelist["fileidtxt"] = str(fileidtxt)
namelist["mergedist"] = str(mergedist)
namelist["reftime"] = str(reftime)
namelist["datadirout"] = str(datadirout)
namelist["fileidout1"] = str(fileidout1)

# Write namelist dictionary to netcdf file for reading 
#  during parallel loop
print("Writing namelist to netcdf file")
nlfileout = Dataset("namelist_po.nc","w",format="NETCDF4")
for k,v in namelist.items():
  setattr(nlfileout, k,  v)
nlfileout.close()

#============================================================
# Get universal information
#============================================================

# Get all original files used for tracking
print("Generating list of files")
datasetFin = Dataset(datadirFin+fileidFin+"00000.nc")
datestart = datasetFin.datestart
dateend = datasetFin.dateend

# Generate a list of filenames with dates to search for
start = dt.datetime.strptime(datestart, "%Y%m%d")
end = dt.datetime.strptime(dateend, "%Y%m%d")
datearr = (start + dt.timedelta(days=x) \
          for x in range(0, (end-start).days))
filen = []
for dateobj in datearr:
  filen.append(datadirIM+fileidIM+ \
  dateobj.strftime("%Y%m%d")+"*")

# Reads directory and filenames
n = 0
for f in filen:
  if n==0:
    filenames = glob.glob(f)
  else:
    filenames = filenames+glob.glob(f)
  n=n+1
lenIdir    = len(datadirIM)

#============================================================
# Read FiT text file
#============================================================

# Read text file
f       = open(datadirFi+fileidFi1+fileidtxt,'r') 
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
if subsetsztm:
  obinds = [list(np.where(FiTids==key)[0]) 
           for key in objs]
  obnt  = list(obids.values())
  obmsz = [max([FiTszs[i] for i in inds]) 
                          for inds in obinds]
  obinds1 = np.intersect1d(
             np.where(np.array(obmsz)>=nsz)[0],
             np.where(np.array(obnt)>=nt)[0])
  objs = [objs[i] for i in obinds1]

# Subset by range of object ids
if subsetobs:
  objs = list(np.intersect1d(
               np.array(np.arange(ob1,ob2)),
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

if serialorparallel==1:
  print("Begin serial loop over objects")
  for o in range(len(objs)):
    dp.driver_processFiTobs(o)

# Parrallel loop over PFs
if serialorparallel==2:
  print("Begin parrallel loop over objects")
  Parallel(n_jobs=njobs)(delayed(dp.driver_processFiTobs)(o) \
    for o in range(len(objs)))

#============================================================
# Final clean up tasks
#============================================================

# Remove namelist and filename files
print("Cleaning up")
os.system("rm filenames_po.txt")
os.system("rm namelist_po.nc")
os.system("rm passinfo_po.nc")
os.system("rm -rf __pycache__")

# End timer
endnow = tm.time()
print("Program took "+str(endnow - startnow)+"s")

#============================================================
# End code
#============================================================
