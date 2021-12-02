import numpy as np
from netCDF4 import Dataset
import warnings
warnings.filterwarnings("ignore")
from scipy.ndimage import uniform_filter as uf
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema as lmax
import datetime as dt
from datetime import datetime
from glob import glob
import sys
import time as tm
from joblib import Parallel, delayed
import os

# Import namelist
from namelist_TIMPS import general as gnl
from namelist_TIMPS import addvars as anl

# Import custom libraries
sys.path.insert(0,gnl.fnsdir)
import time_functions as tfns
import misc_functions as mfns

#==================================================================
# Initialize timer
#==================================================================

startnow = tm.time()

#==================================================================
# Generate list of TIMPS files to process
#==================================================================

# Reads directory and file names
if anl.ssdat:

  print("Subsetting by date")

  # Generate a list of filenames with dates to search for
  print("Generating filenames to search for")
  start = dt.datetime.strptime(anl.date1,"%Y%m%d")
  end = dt.datetime.strptime(anl.date2,"%Y%m%d")
  datearr = (start + dt.timedelta(days=x) \
   for x in range(0,(end-start).days))
  filen = []
  for dateobj in datearr:
    filen.append(f'{gnl.datadirTIMPS}{dateobj.strftime("%m")}/\
{gnl.fileidTIMPS}*_{dateobj.strftime("%Y%m%d")}*.nc')

  # Reads directory and filenames
  print("Getting directory and filenames")
  filenamesrun = []
  for f in filen:
    filenamesrun = filenamesrun+sorted(glob(f))

elif anl.ssobs:

  print("Subsetting by object ids")

  print("Generating filenames to search for")
  obids = [i for i in range(int(anl.obid1),int(anl.obid2))]
  filen = []
  for j in range(1,13):
    for i in range(len(obids)): 
      filen.append(f'{gnl.datadirTIMPS}{str(j).zfill(2)}/\
{gnl.fileidTIMPS}{str(obids[i])}*')

  # Reads directory and filenames
  print("Getting directory and filenames")
  filenamesrun = []
  for f in filen:
    print(f'Adding {str(f)}')
    filenamesrun = filenamesrun+sorted(glob(f))

elif anl.sslst:

  print("Subsetting by files in csv")

  print("Generating list of files")
  filenamesrun = list(np.genfromtxt(anl.lstfn,
   delimiter=',',dtype=str))

else:

  print("No subsetting")

  # Find all files in directory
  print("Generating list of files")
  filenamesrun = []
  for i in range(1,13):
    print(f'{gnl.datadirTIMPS}{str(i).zfill(2)}/{gnl.fileidTIMPS}*')
    filenamesrun = filenamesrun+sorted(glob(
     f'{gnl.datadirTIMPS}{str(i).zfill(2)}/{gnl.fileidTIMPS}*'))

#==================================================================
# Load and process data
#==================================================================

def loadproc_data(inputfile):

  warnings.filterwarnings("ignore")

  # Load data from TIMPS file:
  datasetFin = Dataset(inputfile,'r')
  time = datasetFin.variables["time"][:]
  tunt = datasetFin.variables['time'].units
  area = datasetFin.variables["area"][:]
  vrr  = datasetFin.variables["vrr"][:]

  # Smooth time-series of area and VRR:
  fsize = int(len(area)/5)
  if fsize<5: fsize=5
  areaf = uf(area,size=fsize)
  vrrf = uf(vrr,size=fsize)

  # Get gradient of area and VRR:
  dareaf = np.gradient(areaf,time)
  dvrrf = np.gradient(vrrf,time)

  datasetFin.close()

  return(time,tunt,area,vrr,areaf,vrrf,dareaf,dvrrf)

#==================================================================
# Function to get mature indices from a time-series
#==================================================================

def get_mature_inds(ts,lmxi):
  """
  Get mature stage indices for MCS
  1) Time-series of area or vrr
  2) Indices of local maxima within time-series
  3) Factor of the maxima magnitude for which it can be
      defined as mature stage
    
  Output:
  A list of indices representing the mature stage in the 
   given time-series
  
  """

  # Get magnitude of each maxima
  lmx  = [ts[i] for i in lmxi]
   
  # Get magnitude of minima between each maxima
  sepmin = [min(ts[lmxi[i]:lmxi[i+1]]) 
            for i in range(len(lmxi)-1)]
    
  # Group into sets of maxima
  sep   = [True if ((sepmin[i]<lmx[i]*.75) or 
                    (sepmin[i]<lmx[i+1]*.75)) 
           else False for i in range(len(sepmin))]
  grps = [0]; c=0;
  for i in range(len(sep)): 
    if sep[i]: c=c+1
    grps.extend([c])

  # Select the max in each group
  mxinew = []
  for i in np.unique(grps):
    grpi = np.array(lmxi)[np.where(grps==i)]
    mxinew.extend([grpi[ts[grpi].argmax()]])
  lmx  = [ts[i] for i in mxinew]
    
  # Get minima between each actual maxima
  sepmin = [min(ts[:mxinew[0]])]+\
           [min(ts[mxinew[i]:mxinew[i+1]]) 
            for i in range(len(mxinew)-1)]+\
           [min(ts[mxinew[-1]:])]
  sepmini = [np.where(ts==i)[0].tolist()[0] 
             for i in sepmin]
    
  # Check if each maxima is sufficiently larger than 
  #  the surrounding minima
  mxtf = [True if ((lmx[i]*.75>sepmin[i]) and \
                   (lmx[i]*.75>sepmin[i+1])) 
          else False for i in range(len(lmx))]
  im = [[sepmini[i]+n 
         for n in np.where(ts[
          sepmini[i]:sepmini[i+1]]>lmx[i]*.75)[0].tolist()]
         for i in range(len(mxtf)) if mxtf[i]]

  return([item for sublist in im for item in sublist])

#==================================================================
# Function to check if stages are within certain conditions
#==================================================================

def maxchecks(ts,gmd):
    
  # Ensure growth/decay isn't spuriously defined when a 
  #  maxima is too close to the start of the time-series
  fper = ts[:ts.argmax()+1]
  if min(fper)>max(fper)*.75 or fper[0]>max(fper)*.75:
    for i in np.arange(0,ts.argmax(),1):
      gmd[i]=0
    
  # Ensure growth/decay isn't spuriously defined when 
  # a maxima is too close to the end of the time-series
  lper = ts[ts.argmax():]
  if min(lper)>max(lper)*.75 or lper[-1]>max(lper)*.75:
    for i in np.arange(ts.argmax(),len(ts)-1,1):
      gmd[i]=0

  # Ensure growth or decay doesn't occur until variable 
  #  at least below .75 of a maxima
  for ti in np.arange(ts.argmax(),len(ts),1):
    if ts[ti]>max(ts)*.9:
      if gmd[ti]==1 or gmd[ti]==3: gmd[ti]=0
    else:
      break

  return(gmd)

#==================================================================
# Calculate mature stages
#==================================================================

def get_gmd(areaf,vrrf,dareaf,dvrrf):

  # Calculate all local maxima in filtered VRR and area:
  mxia  = lmax(areaf,np.greater)[0].tolist()
  mxiv  = lmax(vrrf,np.greater)[0].tolist()

  # Instance where there are insufficient time-series to 
  #  calculate the stages
  if len(mxia)<1 or len(mxiv)<1 \
                 or max(areaf)*.5<min(areaf) \
                 or max(vrrf)*.5<min(vrrf):
    print("Time-series is insufficient to calculate MCS phases")
    return(False)

  # Separate multiple local maxima into separate mature stages
  ima = get_mature_inds(areaf,mxia)
  imv = get_mature_inds(vrrf,mxiv)

  # Select only the indices that are mature for both
  im = np.intersect1d(ima,imv)

  # Check against area and vrr change for growth/decay
  gmd = [2 if i in im 
    else 1 if (dareaf[i]>0 and dvrrf[i]>0)
    else 3 if (dareaf[i]<0 and dvrrf[i]<0)
    else 0
    for i in np.arange(0,len(areaf),1)]

  # Checks to ensure growth or decay isn't spuriously 
  #  assigned to a maxima that can't be defined as a maxima
  gmd = maxchecks(areaf,gmd)
  gmd = maxchecks(vrrf,gmd)

  # Ensure each phase lasts longer than just two time-steps
  for i,x in enumerate(gmd):
    if i==0: xnow=x; lenx=1; inds=[0]
    elif i==len(gmd)-1:
      lenx = lenx+1; inds = inds+[i]
      if lenx<4: 
        gmd = [0 if m in inds else n 
               for m,n in enumerate(gmd)]
    elif i>0:
      if xnow==x: lenx = lenx+1; inds = inds+[i]
      else:
        if lenx<4: 
          gmd = [0 if m in inds else n 
                 for m,n in enumerate(gmd)]
        xnow=x; lenx=1; inds=[i]


  # Ensure each growth/decay at least doubles or halves
  for i,x in enumerate(gmd):
    if i==0: xnow=x; inds=[0]
    elif i>0:
      if xnow==x: lenx = lenx+1; inds = inds+[i]
      elif i==len(gmd)-1:
        inds = inds+[i]
        if max(areaf[inds])*.5<min(areaf[inds]):
          gmd = [0 if m in inds else n 
                 for m,n in enumerate(gmd)]
        elif max(vrrf[inds])*.5<min(vrrf[inds]):
          gmd = [0 if m in inds else n 
               for m,n in enumerate(gmd)]
      else:
        if xnow==1 or xnow==3:
          if max(areaf[inds])*.5<min(areaf[inds]):
            gmd = [0 if m in inds else n 
                   for m,n in enumerate(gmd)]
          elif max(vrrf[inds])*.5<min(vrrf[inds]):
            gmd = [0 if m in inds else n 
                   for m,n in enumerate(gmd)]
        xnow=x; lenx=1; inds=[i]

  return(gmd)

#==================================================================
# Run code
#==================================================================

def driver_idstages(fn):

  print(f"Working on {filenamesrun[fn]}")

  # Load required data
  time,tunt,area,vrr,areaf,vrrf,dareaf,dvrrf = loadproc_data(
   filenamesrun[fn])

  # Calculate growth, mature, and decay stages
  gmd1 = get_gmd(areaf,vrrf,dareaf,dvrrf)
  if not gmd1: gmd1 = [0]*len(areaf)

  # Write gmd to file
  fileout = Dataset(filenamesrun[fn],'a')

  description = "An index indicating whether the time is part of \
a growth (gmd=1), mature (gmd=2), or decay (gmd=3) phase for the \
TIMPS. gmd=0 indicates an unidentifiable phase at that time."
  mfns.write_var("gmd","TIMPS Phase",
    description,"time",np.int64,"",fileout,gmd1)

#==================================================================
# Loop over TIMPS files in parrallel
#==================================================================

if anl.serialorparallel==1:
  print("Begin serial loop over objects")
  for fn in range(len(filenamesrun)): driver_idstages(fn)

# Parrallel loop over TIMPS
if anl.serialorparallel==2:
  print("Begin parrallel loop over TIMPS")
  Parallel(n_jobs=anl.njobs)(delayed(driver_idstages)(fn) \
    for fn in range(len(filenamesrun)))

#==================================================================
# Final clean up tasks
#==================================================================

# Remove namelist and filename files
print("Cleaning up")
os.system("rm -rf __pycache__")

# End timer
endnow = tm.time()
print(f'Program took {str(endnow - startnow)}s')


#============================================================
# End
