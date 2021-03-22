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
{gnl.fileidTIMPS}{str(obids[i])}*')

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
  if anl.addmaxrr: maxrainrate  = [np.nan]*len(datakeys)
  if anl.addmeanrr: meanrainrate = [np.nan]*len(datakeys)
  if anl.addmedianrr: medianrainrate = [np.nan]*len(datakeys)
  if anl.addstddevrr: stddevrainrate = [np.nan]*len(datakeys)
  if anl.addskewrr: skewrainrate = [np.nan]*len(datakeys)

  # Pieces
  if anl.addpieces: pieces = [18446744073709551614]*len(datakeys)
  if anl.addpiecesc: 
   pieces_c = [18446744073709551614]*len(datakeys)

  # Area and VRR
  if anl.addarea: area = [np.nan]*len(datakeys)
  if anl.addvrr: volrainrate = [np.nan]*len(datakeys)

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

  # Convective variables
  if anl.addconvrain: is_conv_rain = {}
  if anl.addconvarea: 
    convarea = [np.nan]*len(datakeys)
  if anl.addconvvrr: convvrr = [np.nan]*len(datakeys)

  # Shape variables
  if anl.addaxesshape:
    ellipticity = [np.nan]*len(datakeys)
    axlen       = np.zeros((len(datakeys),2))
    axlen[:,:]  = np.nan
    axang       = np.zeros((len(datakeys),2))
    axang[:,:]  = np.nan
    goodness    = [np.nan]*len(datakeys)
  if anl.addaxesshapec or anl.addpiecesc:
    center_c      = np.zeros((len(datakeys),2))
    center_c[:,:] = np.nan
  if anl.addaxesshapec:
    ellipticity_c = [np.nan]*len(datakeys)
    axlen_c       = np.zeros((len(datakeys),2))
    axlen_c[:,:]  = np.nan
    axang_c       = np.zeros((len(datakeys),2))
    axang_c[:,:]  = np.nan
    goodness_c    = [np.nan]*len(datakeys)
  if anl.addasymmetry: 
    asymmetry_lp = [np.nan]*len(datakeys)
  if anl.addasymmetryc: 
    asymmetry_lp_c = [np.nan]*len(datakeys)
  if anl.addfragmentation or anl.addaxesshape:
    fragmentation  = [np.nan]*len(datakeys)
    solidity       = [np.nan]*len(datakeys)
    connectivity   = [np.nan]*len(datakeys)
  if anl.addfragmentationc or anl.addaxesshapec:
    fragmentation_c  = [np.nan]*len(datakeys)
    solidity_c       = [np.nan]*len(datakeys)
    connectivity_c   = [np.nan]*len(datakeys)
  if anl.adddispersion or anl.addaxesshape:
    dispersion = [np.nan]*len(datakeys)
  if anl.adddispersionc or anl.addaxesshapec:
    dispersion_c = [np.nan]*len(datakeys)

  # Perimeter
  if anl.addperimeter: perimeter_lp = [np.nan]*len(datakeys)

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

#==================================================================
# Calculate maximum rain rate
#==================================================================

    if anl.addmaxrr:

      # Calculate maximum rain rate in PF
      maxrainrate[c] = np.amax(instrainnzk)

#==================================================================
# Calculate mean rain rate
#==================================================================

    if anl.addmeanrr:

      # Calculate mean rain rates but exclude pixels with 
      #  no rain
      meanrainrate[c]   = np.mean(instrainnzk)

#==================================================================
# Calculate median rain rate
#==================================================================

    if anl.addmedianrr:

      # Calculate median rain rates but exclude pixels with 
      #  no rain
      medianrainrate[c] = np.median(instrainnzk)

#==================================================================
# Calculate standard deviation of rain rate
#==================================================================

    if anl.addstddevrr:
      # Calculate standard deviation of rain rates but 
      #  exclude pixels with no rain
      stddevrainrate[c] = np.std(instrainnzk)

#==================================================================
# Calculate standard deviation of rain rate
#==================================================================

    if anl.addskewrr:
      # Calculate standard deviation of rain rates but 
      #  exclude pixels with no rain
      skewrainrate[c] = stats.skew(instrainnzk)

#==================================================================
# Add area 
#==================================================================

    if anl.addarea and not anl.addvrr:

      # Calculate area
      ar = efns.calc_area(lonsnzk,latsnzk,
        float(anl.dx),float(anl.dy))

      # Convert to units of km**2 
      area[c] = ar/1000000

#==================================================================
# Add area and volumetric rain rate information
#==================================================================

    if anl.addarea and anl.addvrr:

      # Calculate area and volumetric rain rate
      ar,vrr = efns.calc_area_and_volrainrate(
        lonsnzk,latsnzk,instrainnzk,float(anl.dx),float(anl.dy))

      # Convert units of km**2 and mm hr**-1 km**2
      area[c] = ar/1000000
      volrainrate[c] = vrr/1000000

#==================================================================
# Add only volumetric rain rate information
#==================================================================

    if anl.addvrr and not anl.addarea:

      # Calculate volumetric rain rate
      vrr = efns.calc_area_and_volrainrate(
        lonsnzk,latsnzk,instrainnzk,
        float(anl.dx),float(anl.dy))[1]

      # Convert to units of mm hr**-1 km**2
      volrainrate[c] = vrr/(1000*2)

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
# Add convective rain flag
#==================================================================

    if anl.addconvrain:

      is_conv_rain[k] = np.where(instrain[k]>anl.convrainthold,1,0)

#==================================================================
# Add convective rain area
#==================================================================
  
    if anl.addconvarea and not anl.addconvvrr:
      # Calculate area
      if len(instrain[k][instrain[k]>anl.convrainthold])>0:
        arc = efns.calc_area(
                    lons[k][instrain[k]>anl.convrainthold],
                    lats[k][instrain[k]>anl.convrainthold],
                    float(anl.dx),float(anl.dy))
      else: arc=0

      # Convert units of km**2 and mm hr**-1 km**2
      convarea[c] = arc/1000000

#==================================================================
# Add convective rain volumetric rain rate
#==================================================================
  
    if anl.addconvvrr and not anl.addconvarea:
      # Calculate area and volumetric rain rate
      if len(instrain[k][instrain[k]>anl.convrainthold])>0:
        vrrc = efns.calc_area_and_volrainrate(
                    lons[k][instrain[k]>anl.convrainthold],
                    lats[k][instrain[k]>anl.convrainthold],
                    instrain[k][instrain[k]>anl.convrainthold],
                    float(anl.dx),float(anl.dy))[1]
      else: vrrc=0

      # Convert units of km**2 and mm hr**-1 km**2
      convvrr[c] = vrrc/1000000

#==================================================================
# Add convective rain area and volumetric rain rate
#==================================================================

    if anl.addconvarea and anl.addconvvrr:
      # Calculate area and volumetric rain rate
      if len(instrain[k][instrain[k]>anl.convrainthold])>0:
        arc,vrrc = efns.calc_area_and_volrainrate(
                    lons[k][instrain[k]>anl.convrainthold],
                    lats[k][instrain[k]>anl.convrainthold],
                    instrain[k][instrain[k]>anl.convrainthold],
                    float(anl.dx),float(anl.dy))
      else:
        arc=0
        vrrc=0

      # Convert units of km**2 and mm hr**-1 km**2
      convarea[c] = arc/1000000
      convvrr[c] = vrrc/1000000

#==================================================================
# Start shape code
#==================================================================

    if anl.addaxesshape or \
       anl.addpieces:

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

     #       # Plotting for test
     #       pfns.plot_pixel_axes_earth(xg,yg,df.values,
     #        center=[dataclon[c],dataclat[c]],
     #        axdir=axang[c,:],axlen=axlen[c,:])
     #else:
     #   print("Basic fit: bad")
     #   # Plotting for test
     #   pfns.plot_pixel_axes_earth(xg,yg,df.values)

#==================================================================
# Convective shape prepwork
#==================================================================

    if anl.addaxesshapec or \
       anl.addpiecesc:

      # Find convective locations
      lonsc = lons[k][instrain[k]>anl.convrainthold]
      latsc = lats[k][instrain[k]>anl.convrainthold]

      # Only if enough convective pixels
      if len(latsc)>float(anl.minshapesizec):

        # Define center of convection
        center_c[c,:] = [efns.periodic_cmass(lonsc),
         np.nanmean(latsc)]

        # assign number of pieces to pieces
        if anl.addpiecesc:
          pieces_c[c] = numLc
          
#==================================================================
# Add major/minor axes shape
#==================================================================

        if anl.addaxesshapec:
            
          # Calculate axes and goodness of fit
          axangc1,axlenc1,goodnessc = \
           sfns.fit_ellipse_svd_earth(lonsc,latsc,center_c[c,:],
           goodness=True,dx=anl.dx,dy=anl.dy)

          # Only continue if fit is good
          if goodnessc>anl.minshapegoodc:

            # Assign data
            goodness_c[c] = goodnessc
            axang_c[c,:] = axangc1
            axlen_c[c,:] = axlenc1
            ellipticity_c[c] = 1-(axlen_c[c,1]/axlen_c[c,0])

          #  # Plotting for test
          #  pfns.plot_pixel_axes_earth(xg,yg,df.values,
          #   center=center_c[c,:],
          #   axdir=axang_c[c,:],axlen=axlen_c[c,:])
          #
          #else:
          #  print("Convective fit: bad ")
          #
          #  # Plotting for test
          #  pfns.plot_pixel_axes_earth(xg,yg,df.values)

      #else: 
      #  print("Not enough convective pixels for a fit")

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
    mfns.write_var("localsolartime","Local solar time",
      description,"time",np.int64,"",fileout,localsolartime)

#==================================================================
# Write max rain rate to file
#==================================================================

  if anl.addmaxrr:

    description = "Maximum rain rate within PF"
    mfns.write_var("maxrainrate","Max rain rate",description,
     "time",np.float64,"mm/hr",fileout,maxrainrate)

#==================================================================
# Write mean rain rate to file
#==================================================================

  if anl.addmeanrr:

    description = "Mean rain rate within PF excluding pixels with zero rain rate"
    mfns.write_var("meanrainrate","Mean rain rate",description,
     "time",np.float64,"mm/hr",fileout,meanrainrate)

#==================================================================
# Write median rain rate to file
#==================================================================

  if anl.addmedianrr:

    description = "Median rain rate within PF excluding pixels with zero rain rate"
    mfns.write_var("medianrainrate","Median rain rate",description,
     "time",np.float64,"mm/hr",fileout,medianrainrate)

#==================================================================
# Write standard deviation rate to file
#==================================================================

  if anl.addstddevrr:

    description = "Standard deviation of rain rates within PF excluding pixels with zero rain rate"
    mfns.write_var("stddevrainrate",
     "Standard deviation of rain rates",description,"time",
      np.float64,"mm/hr",fileout,stddevrainrate)

#==================================================================
# Write standard deviation rate to file
#==================================================================

  if anl.addskewrr:

    description = "Skewness of rain rates within PF excluding pixels with zero rain rate"
    mfns.write_var("skewrainrate",
     "Skewness of rain rates",description,"time",
      np.float64,"",fileout,skewrainrate)

#==================================================================
# Write pieces to file
#==================================================================

  if anl.addpieces:

    description = "Number of disconnected pieces making up the precipitation system"
    mfns.write_var("pieces","Pieces",description,"time",'u8',
     "",fileout,pieces)

#==================================================================
# Write convective pieces to file
#==================================================================

  if anl.addpiecesc:

    description = "Number of disconnected pieces making up the convective component of the system"
    mfns.write_var("pieces_c","Convective pieces",description,
     "time",'u8',"",fileout,pieces_c)

#==================================================================
# Write area to file
#==================================================================

  if anl.addarea:

    description = "Area within PF excluding pixels with zero rain rate"
    mfns.write_var("area","Area",description,"time",np.float64,
     "km^2",fileout,area)

#==================================================================
# Write volumetric rain rate to file
#==================================================================

  if anl.addvrr:

    description = "Volumetric rain rate within PF excluding pixels with zero rain rate"
    mfns.write_var("volrainrate","Volumetric rain rate",
     description,"time",np.float64,"mm hr^-1 km^2",fileout,
     volrainrate)

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
# Write convective information to file
#==================================================================

  if anl.addconvrain or anl.addconvarea or \
    anl.addconvvrr:
    fileout.conv_rain_threshold = anl.convrainthold
 
  if anl.addconvarea:
    description = "Area of locations with rain rates greater than convective rain rate threshold"
    mfns.write_var("area_c","Convective area",
     description,"time",np.float64,"",fileout,convarea)

  if anl.addconvvrr:
    description = "Volumetric rain rate of locations with rain rates greater than convective rain rate threshold"
    mfns.write_var("vrr_c","Convective volumetric rain rate",
    description,"time",np.float64,"",fileout,convvrr)

  if anl.addconvrain:
    format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."
    description = "Binary indicating if a PF location has a convective rain rate (1 = convective, 0 = not convective)"
    mfns.write_group("is_conv_rain",
     "Location has a convective rain rate",description,"",
     format1,fileout,is_conv_rain)

#==================================================================
# Write ellipticity to file
#==================================================================

  if anl.addaxesshape:

    description = "Ellipticity factor. Calculated as 1-(major axis length/minor axis length). 1 = highly elliptical. 0 = spherical."
    mfns.write_var("ellipticity","Ellipticity factor",
      description,"time",np.float64,"",fileout,ellipticity)

#==================================================================
# Write length of axes
#==================================================================

    description = "Length of major axis"
    mfns.write_var("mjrax_length","Major axis length",
     description,"time",np.float64,"m",fileout,axlen[:,0])

    description = "Length of minor axis"
    mfns.write_var("mnrax_length","Minor axis length",
     description,"time",np.float64,"m",fileout,axlen[:,1])

#==================================================================
# Write angle major axis makes from north
#==================================================================

    description = "Angle major axis makes with northward vector"
    mfns.write_var("mjrax_angle","Major axis angle",
     description,"time",np.float64,"degrees",fileout,axang[:,0])

    description = "Angle minor axis makes with northward vector"
    mfns.write_var("mnrax_angle","Minor axis angle",
     description,"time",np.float64,"degrees",fileout,axang[:,1])

#==================================================================
# Write goodness of axes fit for convection
#==================================================================

    description = "Goodness of axes fit to all pixels"
    mfns.write_var("goodness","Axes goodness",
     description,"time",np.float64,"",fileout,goodness)

#==================================================================
# Write ellipticity for convection
#==================================================================

  if anl.addaxesshapec:

    description = "Ellipticity factor for convective pixels. Calculated as 1-(major axis length/minor axis length). 1 = highly elliptical. 0 = spherical."
    mfns.write_var("ellipticity_c","Convective ellipticity factor",
      description,"time",np.float64,"",fileout,ellipticity_c)

#==================================================================
# Write length of axes for convection
#==================================================================

    description = "Length of major axis for convective pixels"
    mfns.write_var("mjrax_length_c","Convective major axis length",
     description,"time",np.float64,"m",fileout,axlen_c[:,0])

    description = "Length of minor axis for convective pixels"
    mfns.write_var("mnrax_length_c","Convective minor axis length",
     description,"time",np.float64,"m",fileout,axlen_c[:,1])

#==================================================================
# Write angle major axis makes from north for convection
#==================================================================

    description = "Angle major axis for convective pixels makes with northward vector"
    mfns.write_var("mjrax_angle_c","Convective major axis angle",
     description,"time",np.float64,"degrees",fileout,axang_c[:,0])

    description = "Angle minor axis for convective pixels makes with northward vector"
    mfns.write_var("mnrax_angle_c","Convective minor axis angle",
     description,"time",np.float64,"degrees",fileout,axang_c[:,1])

#==================================================================
# Write goodness of axes fit for convection
#==================================================================

    description = "Goodness of axes fit to convective pixels"
    mfns.write_var("goodness_c","Convective axes goodness",
     description,"time",np.float64,"",fileout,goodness_c)

#==================================================================
# Write land information to file
#==================================================================

  if anl.addlandinfo:
    description = "1 if center of PF is over land. 0 if not."
    mfns.write_var("cPF_over_land","Center of PF over land",
      description,"time",'u8',"",fileout,cPF_over_land)

    format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."
    description = "Binary indicating if a PF location is over land (1 = over land, 0 = not over land)"
    mfns.write_group("lPF_over_land",
      "Location over land",description,"",
      format1,fileout,lPF_over_land)

#==================================================================
# Write TC information to file
#==================================================================

  if anl.addTCinfo:
    
    if writeTCdata: 

    # If any part of PF IS within TC radius set TC attribute as True and write
    #  other variables describing which parts of PF are within TC radius 
      fileout.within_TC = "True"

      # Write cPF_in_TC
      description = "If PF center is within a maximum radius of TC center this is 1, if not, this is 0."
      mfns.write_var("cPF_in_TC",
        "Proximity of PF center to TC",description,"time",
        'u8',"",fileout,cPF_in_TC)

      # Write dist_cPF_cTC
      description = "Calculated as the geodesic distance of the PF center to TC center"
      mfns.write_var("dist_cPF_cTC",
        "Distance of PF center to TC center",description,
        "time",np.float64,"km",fileout,dist_cPF_cTC)

      # Write TCradius
      description = "Calculated by finding the largest radius of a TC from IBtracs and doubling."
      mfns.write_var("TCrad_cPF","Radius of TC",
        description,"time",np.float64,"km",fileout,TCrad_cPF)

      # Write list of TCname_cPF as an attribute of PF file
      fileout.TCname_cPF = TCname_cPF
     
      # Create or open a group for lPF_in_TC
      format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."
      description = "Binary indicating if a PF location is within a TC (1 = within TCradius, 0 = not within TC radius)"
      mfns.write_group("lPF_in_TC","Location within TC",
        description,"",format1,fileout,lPF_in_TC)

      # Create or open a group for dist_lPF_cTC
      description = "Calculated as geodesic distance to TC center from PF location"
      mfns.write_group("dist_lPF_cTC",
        "Distance to TC center from PF location",description,
        "km",format1,fileout,dist_lPF_cTC)

      # Create or open a group for TCname_lPF
      description = "Name of TC if given location is within TC"
      mfns.write_group("TCname_lPF",
        "Name of TC if given location is within TC",
        description,"",format1,fileout,TCname_lPF)

      # Create or open a group for TCrad_lPF
      description = "Calculated by finding the largest radius of a TC from IBtracs and doubling."
      mfns.write_group("TCrad_lPF","Radius of TC",
        description,"",format1,fileout,TCrad_lPF)

    else:
    # If not, just set TC attribute as False
      fileout.within_TC = "False"

#==================================================================
# Close current file
#==================================================================

  fileout.close()

#==================================================================
# End processing of current PF
#==================================================================

#==================================================================
# Loop over IPF files in parrallel
#==================================================================

if anl.serialorparallel==1:
  print("Begin serial loop over objects")
  for fn in range(len(filenamesrun)): driver_addvars(fn)

# Parrallel loop over PFs
if anl.serialorparallel==2:
  print("Begin parrallel loop over IPFs")
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
