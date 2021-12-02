#==================================================================
# Add variables to IMERG PF netcdf files. Requires 
#  driver_addvars.py be in same directory.
#
# James Russell 2021
#==================================================================

# Import python libraries (do not change)
from netCDF4 import Dataset
import glob
import datetime as dt
import time as tm
from joblib import Parallel, delayed
import os
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
import sys
import math
import scipy.stats as stats
import warnings
warnings.filterwarnings("ignore") 

# Import namelist
from namelist_TIMPS import general as gnl
from namelist_TIMPS import addvars as anl

# Import custom libraries
sys.path.insert(0,gnl.fnsdir)
import earth_functions as efns
import time_functions as tfns
import shape_functions as sfns
import misc_functions as mfns
import plot_functions as pfns

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
    filenamesrun = filenamesrun+sorted(glob.glob(f))

elif anl.ssobs:

  print("Subsetting by object ids")

  print("Generating filenames to search for")
  obids = [i for i in range(int(anl.obid1),int(anl.obid2))]
  filen = []
  for j in range(1,13):
    for i in range(len(obids)): 
      filen.append(f'{gnl.datadirTIMPS}{str(j).zfill(2)}/\
{gnl.fileidTIMPS}{str(obids[i]).zfill(7)}*')
 
  # Reads directory and filenames
  print("Getting directory and filenames")
  filenamesrun = []
  for f in filen:
    print(f'Adding {str(f)}')
    filenamesrun = filenamesrun+sorted(glob.glob(f))

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
    filenamesrun = filenamesrun+sorted(glob.glob(
     f'{gnl.datadirTIMPS}{str(i).zfill(2)}/{gnl.fileidTIMPS}*'))

#==================================================================
# Define function
#==================================================================

def driver_addvars(fn):
  """
  Driver to calculate various variables for a single PF.
  fn corresponds to the line number in filenames_av.txt.
  """

  warnings.filterwarnings("ignore") 

#==================================================================
# Begin script
#==================================================================

  # Get filename for PF
  f = filenamesrun[fn]
  print(f'Working with file {str(f)}')

  # Open file and assign data
  fd = Dataset(f)
  datalat = fd.groups["lats"].groups["data"]
  datalon = fd.groups["lons"].groups["data"]
  datarain = fd.groups["instrain"].groups["data"]
  dataclat = fd.variables["centrallat"][:]
  dataclon = fd.variables["centrallon"][:]
  datacwlat = fd.variables["centlatwgt"][:]
  datacwlon = fd.variables["centlonwgt"][:]
  datadtim = fd.variables["datetime"][:]
  datatim  = fd.variables["time"][:]

  # Get all times 
  datakeys = datalat.__dict__.keys()

#==================================================================
# Preallocate arrays
#==================================================================

  # Rain rates
  if anl.addmaxrr: maxrainrate  = [0]*len(datakeys)
  if anl.addmeanrr: meanrainrate = [0]*len(datakeys)
  if anl.addmedianrr: medianrainrate = [0]*len(datakeys)
  if anl.addstddevrr: stddevrainrate = [np.nan]*len(datakeys)
  if anl.addskewrr: skewrainrate = [np.nan]*len(datakeys)

  if anl.addmeanrr1mmhr: 
   meanrainrate1mmhr = [0]*len(datakeys)
  if anl.addmedianrr1mmhr: 
   medianrainrate1mmhr = [0]*len(datakeys)
  if anl.addstddevrr1mmhr: 
   stddevrainrate1mmhr = [np.nan]*len(datakeys)
  if anl.addskewrr1mmhr: 
   skewrainrate1mmhr = [np.nan]*len(datakeys)

  if anl.addmeanrr10mmhr: 
   meanrainrate10mmhr = [0]*len(datakeys)
  if anl.addmedianrr10mmhr: 
   medianrainrate10mmhr = [0]*len(datakeys)
  if anl.addstddevrr10mmhr: 
   stddevrainrate10mmhr = [np.nan]*len(datakeys)
  if anl.addskewrr10mmhr: 
   skewrainrate10mmhr = [np.nan]*len(datakeys)

  # Pieces
  if anl.addpieces: 
   pieces = [18446744073709551614]*len(datakeys)
  if anl.addpieces1mmhr: 
   pieces1mmhr = [18446744073709551614]*len(datakeys)
  if anl.addpieces10mmhr: 
   pieces10mmhr = [18446744073709551614]*len(datakeys)

  # Area and VRR
  if anl.addarea: area = [0]*len(datakeys)
  if anl.addvrr: volrainrate = [0]*len(datakeys)

  if anl.addarea1mmhr: area1mmhr = [0]*len(datakeys)
  if anl.addvrr1mmhr: volrainrate1mmhr = [0]*len(datakeys)

  if anl.addarea10mmhr: area10mmhr = [0]*len(datakeys)
  if anl.addvrr10mmhr: volrainrate10mmhr = [0]*len(datakeys)

  # Propagation
  if anl.addpropagation:
    propspd = [np.nan]*len(datakeys)
    propdir = [np.nan]*len(datakeys)
    propspdw = [np.nan]*len(datakeys)
    propdirw = [np.nan]*len(datakeys)

  # Time variable
  if anl.addlocaltime: localsolartime = [""]*len(datakeys)

  # TC variables
  if anl.addTCinfo:
    dist_cPF_cTC = [np.nan]*len(datakeys)
    cPF_in_TC    = [18446744073709551614]*len(datakeys)
    TCname_cPF   = [""]*len(datakeys) 
    TCrad_cPF    = [np.nan]*len(datakeys)
    lPF_in_TC    = {}
    dist_lPF_cTC = {}
    TCname_lPF   = {}
    TCrad_lPF    = {}
    writeTCdata   = False

  # Land variable
  if anl.addlandinfo:
    cPF_over_land = [18446744073709551614]*len(datakeys)
    lPF_over_land = {}

  # Shape variables
  if anl.addaxesshape:
    ellipticity = [np.nan]*len(datakeys)
    axlen       = np.zeros((len(datakeys),2))
    axlen[:,:]  = np.nan
    axang       = np.zeros((len(datakeys),2))
    axang[:,:]  = np.nan
    goodness    = [np.nan]*len(datakeys)

  if anl.addaxesshape1mmhr:
    ellipticity_1mmhr = [np.nan]*len(datakeys)
    axlen_1mmhr       = np.zeros((len(datakeys),2))
    axlen_1mmhr[:,:]  = np.nan
    axang_1mmhr       = np.zeros((len(datakeys),2))
    axang_1mmhr[:,:]  = np.nan
    goodness_1mmhr    = [np.nan]*len(datakeys)

  if anl.addaxesshape10mmhr:
    ellipticity_10mmhr = [np.nan]*len(datakeys)
    axlen_10mmhr       = np.zeros((len(datakeys),2))
    axlen_10mmhr[:,:]  = np.nan
    axang_10mmhr       = np.zeros((len(datakeys),2))
    axang_10mmhr[:,:]  = np.nan
    goodness_10mmhr    = [np.nan]*len(datakeys)

#==================================================================
# Begin loop over times
#==================================================================

  c = 0
  for k in datakeys:

    # If first time read, create dictionary
    if c==0:
      lats     = {k:datalat.getncattr(k)}
      lons     = {k:datalon.getncattr(k)}
      instrain = {k:datarain.getncattr(k)}

    # All other times, read data into dictionary
    else:
      lats[k] = datalat.getncattr(k)
      lons[k] = datalon.getncattr(k)
      instrain[k] = datarain.getncattr(k)

    # Get locations of all pixels with nonzero precipitation
    latsnzk = lats[k][instrain[k]>0]
    lonsnzk = lons[k][instrain[k]>0]
    instrainnzk = instrain[k][instrain[k]>0]

    # Get locations of all pixels with nonzero precipitation
    lats1mmhr = lats[k][instrain[k]>1]
    lons1mmhr = lons[k][instrain[k]>1]
    instrain1mmhr = instrain[k][instrain[k]>1]

    # Get locations of all pixels with nonzero precipitation
    lats10mmhr = lats[k][instrain[k]>10]
    lons10mmhr = lons[k][instrain[k]>10]
    instrain10mmhr = instrain[k][instrain[k]>10]

#==================================================================
# Calculate maximum rain rate
#==================================================================

    if anl.addmaxrr:
      maxrainrate[c] = np.amax(instrainnzk)

#==================================================================
# Calculate mean rain rate
#==================================================================

    # Calculate mean rain rates but exclude pixels with no rain
    if anl.addmeanrr and len(instrainnzk)>0:
      meanrainrate[c]   = np.mean(instrainnzk)

    # Calculate mean rain rates but exclude pixels with <1mm/hr
    if anl.addmeanrr1mmhr and len(instrain1mmhr)>0:
      meanrainrate1mmhr[c]   = np.mean(instrain1mmhr)

    # Calculate mean rain rates but exclude pixels with <10mm/hr
    if anl.addmeanrr10mmhr and len(instrain10mmhr)>0:
      meanrainrate10mmhr[c]   = np.mean(instrain10mmhr)

#==================================================================
# Calculate median rain rate
#==================================================================

    # Calculate median rain rates but exclude pixels with no rain
    if anl.addmedianrr and len(instrainnzk)>0:
      medianrainrate[c] = np.median(instrainnzk)

    # Calculate median rain rates but exclude pixels with <1mm/hr
    if anl.addmedianrr1mmhr and len(instrain1mmhr)>0:
      medianrainrate1mmhr[c] = np.median(instrain1mmhr)

    # Calculate median rain rates but exclude pixels with <1mm/hr
    if anl.addmedianrr10mmhr and len(instrain10mmhr)>0:
      medianrainrate10mmhr[c] = np.median(instrain10mmhr)

#==================================================================
# Calculate standard deviation of rain rate
#==================================================================

    # Calculate standard deviation of rain rates but exclude 
    # pixels with no rain
    if anl.addstddevrr and len(instrainnzk)>0:
      stddevrainrate[c] = np.std(instrainnzk)

    # Calculate standard deviation of rain rates but exclude 
    # pixels with <1mm/hr
    if anl.addstddevrr1mmhr and len(instrain1mmhr)>0:
      stddevrainrate1mmhr[c] = np.std(instrain1mmhr)

    # Calculate standard deviation of rain rates but exclude 
    # pixels with <10mm/hr
    if anl.addstddevrr10mmhr and len(instrain10mmhr)>0:
      stddevrainrate10mmhr[c] = np.std(instrain10mmhr)

#==================================================================
# Calculate standard deviation of rain rate
#==================================================================

    # Calculate standard deviation of rain rates but exclude 
    # pixels with no rain
    if anl.addskewrr and len(instrainnzk)>0:
      skewrainrate[c] = stats.skew(instrainnzk)

    # Calculate standard deviation of rain rates but exclude 
    # pixels with no rain
    if anl.addskewrr1mmhr and len(instrain1mmhr)>0:
      skewrainrate1mmhr[c] = stats.skew(instrain1mmhr)

    # Calculate standard deviation of rain rates but exclude 
    # pixels with no rain
    if anl.addskewrr10mmhr and len(instrain10mmhr)>0:
      skewrainrate10mmhr[c] = stats.skew(instrain10mmhr)

#==================================================================
# Add area 
#==================================================================

    # Calculate area of non zero pixels
    if anl.addarea and not anl.addvrr and \
       len(instrainnzk)>0:
      ar = efns.calc_area(lonsnzk,latsnzk,
        float(anl.dx),float(anl.dy))
      area[c] = ar/1000000

    # Calculate area of >1mm/hr pixels
    if anl.addarea1mmhr and not anl.addvrr1mmhr and \
       len(instrain1mmhr)>0:
      ar1mmhr = efns.calc_area(lons1mmhr,lats1mmhr,
        float(anl.dx),float(anl.dy))
      area1mmhr[c] = ar1mmhr/1000000

    # Calculate area of >1mm/hr pixels
    if anl.addarea10mmhr and not anl.addvrr10mmhr and \
       len(instrain10mmhr)>0:
      ar10mmhr = efns.calc_area(lons10mmhr,lats10mmhr,
        float(anl.dx),float(anl.dy))
      area10mmhr[c] = ar10mmhr/1000000

#==================================================================
# Add area and volumetric rain rate information
#==================================================================

    # Calculate area and volumetric rain rate of non-zero pixels
    if anl.addarea and anl.addvrr and \
       len(instrainnzk)>0:
      ar,vrr = efns.calc_area_and_volrainrate(
        lonsnzk,latsnzk,instrainnzk,float(anl.dx),float(anl.dy))
      area[c] = ar/1000000
      volrainrate[c] = vrr/1000000

    # Calculate area and volumetric rain rate of >1mm/hr pixels
    if anl.addarea1mmhr and anl.addvrr1mmhr and \
       len(instrain1mmhr)>0:
      ar1mmhr,vrr1mmhr = efns.calc_area_and_volrainrate(
        lons1mmhr,lats1mmhr,instrain1mmhr,
        float(anl.dx),float(anl.dy))
      area1mmhr[c] = ar1mmhr/1000000
      volrainrate1mmhr[c] = vrr1mmhr/1000000

    # Calculate area and volumetric rain rate of >1mm/hr pixels
    if anl.addarea10mmhr and anl.addvrr10mmhr and \
       len(instrain10mmhr)>0:
      ar10mmhr,vrr10mmhr = efns.calc_area_and_volrainrate(
        lons10mmhr,lats10mmhr,instrain10mmhr,
        float(anl.dx),float(anl.dy))
      area10mmhr[c] = ar10mmhr/1000000
      volrainrate10mmhr[c] = vrr10mmhr/1000000

#==================================================================
# Add only volumetric rain rate information
#==================================================================

    # Calculate volumetric rain rate of non-zero pixels
    if anl.addvrr and not anl.addarea and \
       len(instrainnzk)>0:
      vrr = efns.calc_area_and_volrainrate(
        lonsnzk,latsnzk,instrainnzk,
        float(anl.dx),float(anl.dy))[1]
      volrainrate[c] = vrr/1000000

    # Calculate volumetric rain rate of >1mm/hr pixels
    if anl.addvrr1mmhr and not anl.addarea1mmhr and \
       len(instrain1mmhr)>0:
      vrr1mmhr = efns.calc_area_and_volrainrate(
        lons1mmhr,lats1mmhr,instrain1mmhr,
        float(anl.dx),float(anl.dy))[1]
      volrainrate1mmhr[c] = vrr1mmhr/1000000

    # Calculate volumetric rain rate of >10mm/hr pixels
    if anl.addvrr10mmhr and not anl.addarea10mmhr and \
       len(instrain10mmhr)>0:
      vrr10mmhr = efns.calc_area_and_volrainrate(
        lons10mmhr,lats10mmhr,instrain10mmhr,
        float(anl.dx),float(anl.dy))[1]
      volrainrate10mmhr[c] = vrr10mmhr/1000000

#==================================================================
# Add local solar time
#==================================================================

    if anl.addlocaltime:
    
      # Calculate local time based on central longitude
      localsolartime[c] = tfns.calc_local_solar_date(
       str(datadtim[c])[0:4],str(datadtim[c])[4:6],
       str(datadtim[c])[6:8],str(datadtim[c])[8:10],
       str(datadtim[c])[10:12],str(0),dataclon[c])

#==================================================================
# Add propagation
#==================================================================

    if anl.addpropagation:
    # Calculate speed and direction of motion of PF centroid

      # Account for objects that aren't more than one time
      if len(dataclat)<2:

        # If only one time in PF, set as 0
        propspd[c] = 0
        propdir[c] = 0
        propspdw[c] = 0
        propdirw[c] = 0

      else:

        # Forward difference for first time
        if c==0:
          propspd[c],propdir[c] = efns.calc_propagation(
            str(datadtim[c]*100),  dataclon[c],  dataclat[c],
            str(datadtim[c+1]*100),dataclon[c+1],dataclat[c+1])

          propspdw[c],propdirw[c] = efns.calc_propagation(
            str(datadtim[c]*100),  datacwlon[c],  datacwlat[c],
            str(datadtim[c+1]*100),datacwlon[c+1],datacwlat[c+1])

        # Backward difference for last time
        elif c==len(dataclon)-1:
          propspd[c],propdir[c] = efns.calc_propagation(
            str(datadtim[c-1]*100),dataclon[c-1],dataclat[c-1],
            str(datadtim[c]*100),  dataclon[c],  dataclat[c])

          propspdw[c],propdirw[c] = efns.calc_propagation(
            str(datadtim[c-1]*100),dataclon[c-1],dataclat[c-1],
            str(datadtim[c]*100),  dataclon[c],  dataclat[c])

        # Centered difference for all other times
        else:
          propspd[c],propdir[c] = efns.calc_propagation(
            str(datadtim[c-1]*100),dataclon[c-1],dataclat[c-1],
            str(datadtim[c+1]*100),dataclon[c+1],dataclat[c+1])

          propspdw[c],propdirw[c] = efns.calc_propagation(
            str(datadtim[c-1]*100),datacwlon[c-1],datacwlat[c-1],
            str(datadtim[c+1]*100),datacwlon[c+1],datacwlat[c+1])

#==================================================================
# Start shape code
#==================================================================

    if (anl.addaxesshape or anl.addpieces) and \
       len(instrainnzk)>0:

      # Generate dataframe for object
      df = mfns.create_2d_dataframe(lonsnzk,latsnzk,
       anl.dx,anl.dy,instrainnzk)
      ny = [float(i) for i in df.index]
      nx = [float(i) for i in df.columns]
      xg, yg = np.meshgrid(nx,ny)

      # assign number of pieces to pieces
      if anl.addpieces:
        numL = sfns.label_wdiags(df)[1]
        pieces[c] = numL
            
#==================================================================
# Add major/minor axes shape
#==================================================================
        
      # Only calculate if there are enough points
      if len(latsnzk)>float(anl.minshapesize):

        # Fragmentation prepwork
        if anl.addaxesshape:
          
          # Calculate axes and goodness
          axang1,axlen1,goodness1 = \
           sfns.fit_ellipse_svd_earth(lonsnzk,latsnzk,
           [dataclon[c],dataclat[c]],goodness=True,
           dx=anl.dx,dy=anl.dy)

          # Only continue if fit is good
          if goodness1>anl.minshapegood:
            
            # Assign data
            goodness[c] = goodness1
            axang[c,:] = axang1
            axlen[c,:] = axlen1
            ellipticity[c] = 1-(axlen[c,1]/axlen[c,0])

#==================================================================
# Start shape code
#==================================================================

    if anl.addaxesshape1mmhr or anl.addpieces1mmhr and \
       len(instrain1mmhr)>0:

      # Generate dataframe for object
      df = mfns.create_2d_dataframe(lons1mmhr,lats1mmhr,
       anl.dx,anl.dy,instrain1mmhr)
      ny = [float(i) for i in df.index]
      nx = [float(i) for i in df.columns]
      xg, yg = np.meshgrid(nx,ny)

      # assign number of pieces to pieces
      if anl.addpieces1mmhr:
        numL = sfns.label_wdiags(df)[1]
        pieces1mmhr[c] = numL
            
#==================================================================
# Add major/minor axes shape
#==================================================================
        
      # Only calculate if there are enough points
      if len(lats1mmhr)>float(anl.minshapesize):

        # Fragmentation prepwork
        if anl.addaxesshape1mmhr:
          
          # Calculate axes and goodness
          axang1,axlen1,goodness1 = \
           sfns.fit_ellipse_svd_earth(lons1mmhr,lats1mmhr,
           [dataclon[c],dataclat[c]],goodness=True,
           dx=anl.dx,dy=anl.dy)

          # Only continue if fit is good
          if goodness1>anl.minshapegood1mmhr:
            
            # Assign data
            goodness_1mmhr[c] = goodness1
            axang_1mmhr[c,:] = axang1
            axlen_1mmhr[c,:] = axlen1
            ellipticity_1mmhr[c] = 1-(
             axlen_1mmhr[c,1]/axlen_1mmhr[c,0])

#==================================================================
# Start shape code
#==================================================================

    if anl.addaxesshape10mmhr or anl.addpieces10mmhr and \
       len(instrain10mmhr)>0:

      # Generate dataframe for object
      df = mfns.create_2d_dataframe(lons10mmhr,lats10mmhr,
       anl.dx,anl.dy,instrain10mmhr)
      ny = [float(i) for i in df.index]
      nx = [float(i) for i in df.columns]
      xg, yg = np.meshgrid(nx,ny)

      # assign number of pieces to pieces
      if anl.addpieces10mmhr:
        numL = sfns.label_wdiags(df)[1]
        pieces10mmhr[c] = numL
            
#==================================================================
# Add major/minor axes shape
#==================================================================
        
      # Only calculate if there are enough points
      if len(lats10mmhr)>float(anl.minshapesize10mmhr):

        # Fragmentation prepwork
        if anl.addaxesshape10mmhr:
          
          # Calculate axes and goodness
          axang1,axlen1,goodness1 = \
           sfns.fit_ellipse_svd_earth(lons10mmhr,lats10mmhr,
           [dataclon[c],dataclat[c]],goodness=True,
           dx=anl.dx,dy=anl.dy)

          # Only continue if fit is good
          if goodness1>anl.minshapegood10mmhr:
            
            # Assign data
            goodness_10mmhr[c] = goodness1
            axang_10mmhr[c,:] = axang1
            axlen_10mmhr[c,:] = axlen1
            ellipticity_10mmhr[c] = 1-(
             axlen_10mmhr[c,1]/axlen_10mmhr[c,0])

#==================================================================
# Add land surface information
#==================================================================

    if anl.addlandinfo:
      
      # Import specific libraries
      import shapely.geometry as sgeom
      
      # Get land shape file
      land = efns.load_land()

      # Check if center of PF over land
      if land.contains(sgeom.Point(dataclon[c], dataclat[c])):
        cPF_over_land[c] = 1
      else:
        cPF_over_land[c] = 0

      # If only a single pixel
      if not hasattr(lons[k], "__len__"):
        if land.contains(sgeom.Point(lons[k],lats[k])):
          latlonland = 1
        else:
          latlonland = 0

      # If larger 
      else:
        latlonland = [0]*len(lons[k])
      
        # Loop over all locations of PF at time k
        for l in range(len(lons[k])):          
      
          # Check if center of PF over land
          if land.contains(sgeom.Point(lons[k][l],lats[k][l])):
            latlonland[l] = 1
          else:
            latlonland[l] = 0

        # Assign scalar/list to key (time) in dictionary
        lPF_over_land[k] = latlonland
     
#==================================================================
# Add TC information
#==================================================================

    if anl.addTCinfo:
      
      # Read tropical cyclone data
      fTC = Dataset(anl.dataTCdir+anl.fileTCid) 

      # Interpolate times and latitudes and get radius
      TCinfo = efns.interp_TC(datadtim[c],fTC)

      if TCinfo["TC"]:

        # Check if center is close to TC
        dist_cPF_cTC[c],cPF_in_TC[c],TCname_cPF[c],TCrad_cPF[c] \
          = efns.calc_if_TC(dataclon[c],dataclat[c],TCinfo,fTC)

        # Ensure variables are not scalar
        if not hasattr(lats[k], "__len__"):
          lats[k] = [lats[k]]
          lons[k] = [lons[k]]

        # Preallocate arrays for TC information to location in PF
        TCdistnow = [0.]*len(lats[k])
        TCbinnow  = [0]*len(lats[k])
        TCnamenow = [""]*len(lats[k])
        TCradnow  = [0.]*len(lats[k])

        # Check all locations in PF for distance to TC center
        for il in range(len(lats[k])):
          TCdistnow[il],TCbinnow[il],TCnamenow[il],TCradnow[il] \
            = efns.calc_if_TC(lons[k][il],lats[k][il],TCinfo,fTC)

        # Assign to dictionaries
        lPF_in_TC[k]    = TCbinnow
        dist_lPF_cTC[k] = TCdistnow
        TCname_lPF[k]   = TCnamenow
        TCrad_lPF[k]    = TCradnow

        # Set boolean to indicate if writing data is required
        if (cPF_in_TC[c]==1) or (1 in lPF_in_TC[k]):
          writeTCdata = True

      else:
        writeTCdata = False

#==================================================================
# End loops over objects and times
#==================================================================
      
    # Advance counter
    c = c + 1

  fd.close()

#==================================================================
# Open file to write data
#==================================================================

  fileout = Dataset(f,'a')

#==================================================================
# Write local time to file
#==================================================================

  if anl.addlocaltime:

    localsolartime = [int(i) for i in localsolartime]
    description = "Calculated as the UTC time plus an offset based on the longitude. The offset is calculated by multiplying the longitude by 24/360. Note: This is not the actual local time. This should typically only be used to calculate times for the diurnal cycle."
    mfns.write_var("LST","Local solar time",
      description,"time",np.int64,"",fileout,localsolartime)

#==================================================================
# Write max rain rate to file
#==================================================================

  if anl.addmaxrr:
    description = "Maximum rain rate of TIMPS"
    mfns.write_var("maxrr","Max rain rate",description,
     "time",np.float64,"mm/hr",fileout,maxrainrate)

#==================================================================
# Write mean rain rate to file
#==================================================================

  if anl.addmeanrr:
    description = "Mean rain rate of non-zero pixels"
    mfns.write_var("meanrr","Mean rain rate",description,
     "time",np.float64,"mm/hr",fileout,meanrainrate)

  if anl.addmeanrr1mmhr:
    description = "Mean rain rate of >1mm/hr pixels"
    mfns.write_var("meanrr1","Mean rain rate (>1mm/hr)",
     description,"time",np.float64,"mm/hr",fileout,
     meanrainrate1mmhr)

  if anl.addmeanrr10mmhr:
    description = "Mean rain rate of >10mm/hr pixels"
    mfns.write_var("meanrr10",
     "Mean rain rate (>10mm/hr)",description,"time",np.float64,
     "mm/hr",fileout,meanrainrate10mmhr)

#==================================================================
# Write median rain rate to file
#==================================================================

  if anl.addmedianrr:
    description = "Median rain rate of non-zero pixels"
    mfns.write_var("medrr","Median rain rate",
     description,"time",np.float64,"mm/hr",fileout,medianrainrate)

  if anl.addmedianrr1mmhr:
    description = "Median rain rate of >1mm/hr pixels"
    mfns.write_var("medrr1",
     "Median rain rate (>1mm/hr)",description,"time",
     np.float64,"mm/hr",fileout,medianrainrate1mmhr)

  if anl.addmedianrr10mmhr:
    description = "Median rain rate of >10mm/hr pixels"
    mfns.write_var("medrr10",
     "Median rain rate (>10mm/hr)",description,"time",
     np.float64,"mm/hr",fileout,medianrainrate10mmhr)

#==================================================================
# Write standard deviation rate to file
#==================================================================

  if anl.addstddevrr:
    description = "Standard deviation of non-zero rain rates"
    mfns.write_var("sdevrr",
     "Standard deviation of rain rates",description,
     "time",np.float64,"mm/hr",fileout,stddevrainrate)

  if anl.addstddevrr1mmhr:
    description = "Standard deviation of >1mm/hr rain rates"
    mfns.write_var("sdevrr1",
     "Standard deviation of rain rates (>1mm/hr)",description,
     "time",np.float64,"mm/hr",fileout,stddevrainrate1mmhr)

  if anl.addstddevrr10mmhr:
    description = "Standard deviation of >10mm/hr rain rates"
    mfns.write_var("sdevrr10",
     "Standard deviation of rain rates (>10mm/hr)",description,
     "time",np.float64,"mm/hr",fileout,stddevrainrate10mmhr)

#==================================================================
# Write standard deviation rate to file
#==================================================================

  if anl.addskewrr:
    description = "Skewness of non-zero rain rates"
    mfns.write_var("skewrr",
     "Skewness of rain rates",description,"time",
      np.float64,"",fileout,skewrainrate)

  if anl.addskewrr1mmhr:
    description = "Skewness of >1mm/hr rain rates"
    mfns.write_var("skewrr1",
     "Skewness of rain rates (>1mm/hr)",description,"time",
      np.float64,"",fileout,skewrainrate1mmhr)

  if anl.addskewrr10mmhr:
    description = "Skewness of >10mm/hr rain rates"
    mfns.write_var("skewrr10",
     "Skewness of rain rates (>10mm/hr)",description,"time",
      np.float64,"",fileout,skewrainrate10mmhr)

#==================================================================
# Write pieces to file
#==================================================================

  if anl.addpieces:

    description = "Number of disconnected pieces making up the TIMPS"
    mfns.write_var("pieces","Pieces",description,"time",'u8',
     "",fileout,pieces)

  if anl.addpieces1mmhr:

    description = "Number of disconnected pieces making up the convective component of the system"
    mfns.write_var("pieces1","Pieces (>1mm/hr)",
     description,"time",'u8',"",fileout,pieces1mmhr)

  if anl.addpieces10mmhr:

    description = "Number of disconnected pieces making up the convective component of the system"
    mfns.write_var("pieces10","Pieces (>10mm/hr)",
     description,"time",'u8',"",fileout,pieces10mmhr)

#==================================================================
# Write area to file
#==================================================================

  if anl.addarea:
    description = "Area of non-zero pixels"
    mfns.write_var("area","Area",description,
     "time",np.float64,"km^2",fileout,area)

  if anl.addarea1mmhr:
    description = "Area of >1mm/hr pixels"
    mfns.write_var("area1","Area (>1mm/hr)",description,
     "time",np.float64,"km^2",fileout,area1mmhr)

  if anl.addarea10mmhr:
    description = "Area of >10mm/hr pixels"
    mfns.write_var("area10","Area (>10mm/hr)",description,
     "time",np.float64,"km^2",fileout,area10mmhr)

#==================================================================
# Write volumetric rain rate to file
#==================================================================

  if anl.addvrr:
    description = "Volumetric rain rate of non-zero pixels"
    mfns.write_var("vrr","Volumetric rain rate",
     description,"time",np.float64,"mm hr^-1 km^2",fileout,
     volrainrate)

  if anl.addvrr1mmhr:
    description = "Volumetric rain rate of >1mm/hr pixels"
    mfns.write_var("vrr1",
     "Volumetric rain rate (>1mmhr)",description,"time",
     np.float64,"mm hr^-1 km^2",fileout,volrainrate1mmhr)

  if anl.addvrr10mmhr:
    description = "Volumetric rain rate of >10mm/hr pixels"
    mfns.write_var("vrr10",
     "Volumetric rain rate (>10mmhr)",description,"time",
     np.float64,"mm hr^-1 km^2",fileout,volrainrate10mmhr)

#==================================================================
# Write propagation to file
#==================================================================

  if anl.addpropagation:

    description = "Calculated as the geodesic distance travelled by centroid divided by time"
    mfns.write_var("propspd","Propagation speed",description,
     "time",np.float64,"m/s",fileout,propspd)

    description = "Calculated as direction centroid is moving toward from north (clockwise)"
    mfns.write_var("propdir","Propagation direction",description,
     "time",np.float64,"degrees",fileout,propdir)

    description = "Calculated as the geodesic distance travelled by precipitation weighted centroid divided by time"
    mfns.write_var("propspdw","Propagation speed (weighted)",
     description,"time",np.float64,"m/s",fileout,propspdw)

    description = "Calculated as direction precipitation weighted centroid is moving toward from north (clockwise)"
    mfns.write_var("propdirw","Propagation direction (weighted)",
     description,"time",np.float64,"degrees",fileout,propdirw)

#==================================================================
# Write axes parameters for >1mm/hr pixels
#==================================================================

  if anl.addaxesshape:

    description = "Ellipticity factor for non-zero pixels. Calculated as 1-(major axis length/minor axis length). 1 = highly elliptical. 0 = spherical."
    mfns.write_var("ell",
     "Ellipticity factor",description,"time",
     np.float64,"",fileout,ellipticity)

    description = "Length of major axis for non-zero pixels"
    mfns.write_var("mjraxlen",
     "Major axis length",description,"time",
     np.float64,"m",fileout,axlen[:,0])

    description = "Length of minor axis for non-zero pixels"
    mfns.write_var("mnraxlen",
     "Minor axis length",description,"time",
     np.float64,"m",fileout,axlen[:,1])

    description = "Angle major axis for non-zero makes with northward vector"
    mfns.write_var("mjraxang",
     "Major axis angle",description,"time",
     np.float64,"degrees",fileout,axang[:,0])

    description = "Angle minor axis for non-zero pixels makes with northward vector"
    mfns.write_var("mnraxangle",
     "Minor axis angle",description,"time",
     np.float64,"degrees",fileout,axang[:,1])

    description = "Goodness of axes fit to non-zero pixels"
    mfns.write_var("goodax","Axes goodness",
     description,"time",np.float64,"",fileout,goodness)

#==================================================================
# Write axes parameters for >1mm/hr pixels
#==================================================================

  if anl.addaxesshape1mmhr:

    description = "Ellipticity factor for >1mm/hr pixels. Calculated as 1-(major axis length/minor axis length). 1 = highly elliptical. 0 = spherical."
    mfns.write_var("ell1",
     "Ellipticity factor (>1mm/hr)",description,"time",
     np.float64,"",fileout,ellipticity1mmhr)

    description = "Length of major axis for >1mm/hr pixels"
    mfns.write_var("mjraxlen1",
     "Major axis length (>1mm/hr)",description,"time",
     np.float64,"m",fileout,axlen1mmhr[:,0])

    description = "Length of minor axis for >1mm/hr pixels"
    mfns.write_var("mnraxlen1",
     "Minor axis length (>1mm/hr)",description,"time",
     np.float64,"m",fileout,axlen1mmhr[:,1])

    description = "Angle major axis for >1mm/hr makes with northward vector"
    mfns.write_var("mjraxang1",
     "Major axis angle (>1mm/hr)",description,"time",
     np.float64,"degrees",fileout,axang1mmhr[:,0])

    description = "Angle minor axis for >1mm/hr pixels makes with northward vector"
    mfns.write_var("mnraxangle1",
     "Minor axis angle (>1mm/hr)",description,"time",
     np.float64,"degrees",fileout,axang1mmhr[:,1])

    description = "Goodness of axes fit to >1mm/hr pixels"
    mfns.write_var("goodax1","Axes goodness (>1mm/hr)",
     description,"time",np.float64,"",fileout,goodness1mmhr)

#==================================================================
# Write axes parameters for >10mm/hr pixels
#==================================================================

  if anl.addaxesshape10mmhr:

    description = "Ellipticity factor for >10mm/hr pixels. Calculated as 1-(major axis length/minor axis length). 1 = highly elliptical. 0 = spherical."
    mfns.write_var("ell10",
     "Ellipticity factor (>10mm/hr)",description,"time",
     np.float64,"",fileout,ellipticity10mmhr)

    description = "Length of major axis for >10mm/hr pixels"
    mfns.write_var("mjraxlen10",
     "Major axis length (>10mm/hr)",description,"time",
     np.float64,"m",fileout,axlen10mmhr[:,0])

    description = "Length of minor axis for >10mm/hr pixels"
    mfns.write_var("mnraxlen10",
     "Minor axis length (>10mm/hr)",description,"time",
     np.float64,"m",fileout,axlen10mmhr[:,1])

    description = "Angle major axis for >10mm/hr makes with northward vector"
    mfns.write_var("mjraxang10",
     "Major axis angle (>10mm/hr)",description,"time",
     np.float64,"degrees",fileout,axang10mmhr[:,0])

    description = "Angle minor axis for >10mm/hr pixels makes with northward vector"
    mfns.write_var("mnraxangle10",
     "Minor axis angle (>10mm/hr)",description,"time",
     np.float64,"degrees",fileout,axang10mmhr[:,1])

    description = "Goodness of axes fit to >10mm/hr pixels"
    mfns.write_var("goodax10","Axes goodness (>10mm/hr)",
     description,"time",np.float64,"",fileout,goodness10mmhr)

#==================================================================
# Write land information to file
#==================================================================

  if anl.addlandinfo:
    description = "1 if center of TIMPS is over land. 0 if not."
    mfns.write_var("cent_over_land","Center of TIMPS over land",
      description,"time",'u8',"",fileout,cPF_over_land)

    format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."
    description = "Binary indicating if a TIMPS location is over land (1 = over land, 0 = not over land)"
    mfns.write_group("loc_over_land",
      "Location over land",description,"",
      format1,fileout,lPF_over_land,f)

#==================================================================
# Write TC information to file
#==================================================================

  if anl.addTCinfo:
    
    if writeTCdata: 

    # If any part of TIMPS IS within TC radius set TC attribute 
    #  as True and write other variables describing which parts of
    #  TIMPS are within TC radius 
      fileout.within_TC = "True"

      # Write cent_in_TC
      description = "If TIMPS center is within a maximum radius of TC center this is 1, if not, this is 0."
      mfns.write_var("cent_in_TC",
        "Proximity of TIMPS center to TC",description,"time",
        'u8',"",fileout,cPF_in_TC)

      # Write dist_cent_cTC
      description = "Calculated as the geodesic distance of the TIMPS center to TC center"
      mfns.write_var("dist_cent_cTC",
        "Distance of TIMPS center to TC center",description,
        "time",np.float64,"km",fileout,dist_cPF_cTC)

      # Write TCradius
      description = "Calculated by finding the largest radius of a TC from IBtracs and doubling."
      mfns.write_var("TCrad_cent","Radius of TC",
        description,"time",np.float64,"km",fileout,TCrad_cPF)

      # Write list of TCname_cent as an attribute of TIMPS file
      fileout.TCname_cent = TCname_cPF
     
      # Create or open a group for loc_in_TC
      format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."
      description = "Binary indicating if a TIMPS location is within a TC (1 = within TCradius, 0 = not within TC radius)"
      mfns.write_group("loc_in_TC","Location within TC",
        description,"",format1,fileout,lPF_in_TC,f)

      # Create or open a group for dist_loc_cTC
      description = "Calculated as geodesic distance to TC center from TIMPS location"
      mfns.write_group("dist_loc_cTC",
        "Distance to TC center from TIMPS location",description,
        "km",format1,fileout,dist_lPF_cTC,f)

      # Create or open a group for TCname_loc
      description = "Name of TC if given location is within TC"
      mfns.write_group("TCname_loc",
        "Name of TC if given location is within TC",
        description,"",format1,fileout,TCname_lPF,f)

      # Create or open a group for TCrad_loc
      description = "Calculated by finding the largest radius of a TC from IBtracs and doubling."
      mfns.write_group("TCrad_loc","Radius of TC",
        description,"",format1,fileout,TCrad_lPF,f)

    else:
    # If not, just set TC attribute as False
      fileout.within_TC = "False"

#==================================================================
# Close current file
#==================================================================

  fileout.close()

#==================================================================
# End processing of current TIMPS
#==================================================================

#==================================================================
# Loop over TIMPS files in parrallel
#==================================================================

if anl.serialorparallel==1:
  print("Begin serial loop over objects")
  for fn in range(len(filenamesrun)): driver_addvars(fn)

# Parrallel loop over TIMPS
if anl.serialorparallel==2:
  print("Begin parrallel loop over TIMPS")
  Parallel(n_jobs=anl.njobs)(delayed(driver_addvars)(fn) \
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

#==================================================================
# End code
#==================================================================
