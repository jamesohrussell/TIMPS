#==================================================================
# Add variables to IMERG PF netcdf files. Requires 
#  driver_addvars.py be in same directory.
#
# James Russell 2019
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
import geophys_functions as gfns
import time_functions as tfns
import shape_functions as sfns
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
    filen.append(f'{gnl.datadirTIPS}{dateobj.strftime("%m")}/\
{gnl.fileidTIPS}*_{dateobj.strftime("%Y%m%d")}*.nc')

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
      filen.append(f'{gnl.datadirTIPS}{str(j).zfill(2)}/\
{gnl.fileidTIPS}{str(obids[i])}*')

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
    filenamesrun = filenamesrun+sorted(glob.glob(
     f'{gnl.datadirTIPS}{str(i).zfill(2)}/{gnl.fileidTIPS}*'))

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
  if anl.addaxesshapec or anl.addpiecesc:
    center_c      = np.zeros((len(datakeys),2))
    center_c[:,:] = np.nan
  if anl.addaxesshapec:
    ellipticity_c = [np.nan]*len(datakeys)
    axlen_c       = np.zeros((len(datakeys),2))
    axlen_c[:,:]  = np.nan
    axang_c       = np.zeros((len(datakeys),2))
    axang_c[:,:]  = np.nan
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
      ar = gfns.calc_area(lonsnzk,latsnzk,
        float(anl.dx),float(anl.dy))

      # Convert to units of km**2 
      area[c] = ar/1000000

#==================================================================
# Add area and volumetric rain rate information
#==================================================================

    if anl.addarea and anl.addvrr:

      # Calculate area and volumetric rain rate
      ar,vrr = gfns.calc_area_and_volrainrate(
        lonsnzk,latsnzk,instrainnzk,float(anl.dx),float(anl.dy))

      # Convert units of km**2 and mm hr**-1 km**2
      area[c] = ar/1000000
      volrainrate[c] = vrr/1000000

#==================================================================
# Add only volumetric rain rate information
#==================================================================

    if anl.addvrr and not anl.addarea:

      # Calculate volumetric rain rate
      vrr = gfns.calc_area_and_volrainrate(
        lonsnzk,latsnzk,instrainnzk,
        float(anl.dx),float(anl.dy))[1]

      # Convert to units of mm hr**-1 km**2
      volrainrate[c] = vrr/(1000*2)

#==================================================================
# Add local solar time
#==================================================================

    if anl.addlocaltime:
    
      # Calculate local time based on central longitude
      localsolartime[c] = tfns.calc_local_solar_time(
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
          propspd[c],propdir[c] = gfns.calc_propagation(
            str(datadtim[c]*100),  dataclon[c],  dataclat[c],
            str(datadtim[c+1]*100),dataclon[c+1],dataclat[c+1])

          propspdw[c],propdirw[c] = gfns.calc_propagation(
            str(datadtim[c]*100),  datacwlon[c],  datacwlat[c],
            str(datadtim[c+1]*100),datacwlon[c+1],datacwlat[c+1])

        # Backward difference for last time
        elif c==len(dataclon)-1:
          propspd[c],propdir[c] = gfns.calc_propagation(
            str(datadtim[c-1]*100),dataclon[c-1],dataclat[c-1],
            str(datadtim[c]*100),  dataclon[c],  dataclat[c])

          propspdw[c],propdirw[c] = gfns.calc_propagation(
            str(datadtim[c-1]*100),dataclon[c-1],dataclat[c-1],
            str(datadtim[c]*100),  dataclon[c],  dataclat[c])

        # Centered difference for all other times
        else:
          propspd[c],propdir[c] = gfns.calc_propagation(
            str(datadtim[c-1]*100),dataclon[c-1],dataclat[c-1],
            str(datadtim[c+1]*100),dataclon[c+1],dataclat[c+1])

          propspdw[c],propdirw[c] = gfns.calc_propagation(
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
        arc = gfns.calc_area(
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
        vrrc = gfns.calc_area_and_volrainrate(
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
        arc,vrrc = gfns.calc_area_and_volrainrate(
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

    if anl.addperimeter or \
       anl.addasymmetry or \
       anl.addfragmentation or \
       anl.addaxesshape or \
       anl.adddispersion or \
       anl.addpieces:

      # Generate dataframe for object
      df = mfns.create_2d_dataframe(lonsnzk,latsnzk,
                                   anl.dx,anl.dy,instrainnzk)

      # Find and label all contiguous areas within object
      labels, numL = sfns.label_wdiags(df)
      ny = [float(i) for i in df.index]
      nx = [float(i) for i in df.columns]
      df1 = pd.DataFrame(labels,[str(i) for i in ny],
                                [str(i) for i in nx])

      # assign number of pieces to pieces
      if anl.addpieces:
        pieces[c] = numL

      # Only calculate if there are enough points
      if len(latsnzk)>float(anl.minshapesize):

        # Fragmentation prepwork
        if anl.addfragmentation or \
           anl.addaxesshape:

          # Predefine solidity array
          solid = [0.]*numL
          
          # Define points in the grid
          xg, yg = np.meshgrid(nx,ny)
          points = [[round(float(x),2) for x in y] 
           for y in np.vstack((xg.flatten(),yg.flatten())).T]

        # Loop over all separate pieces
        areas    = [0.]*numL
        centersx = [0.]*numL
        centersy = [0.]*numL
        for i in range(1,numL+1):
        
          # Line 2: Identify pairs of indices for current piece
          # Line 1: Convert all pairs to floats
          indpairs = [(round(float(z[1]),2),round(float(z[0]),2))
           for z in df1[df1==i].stack().index.tolist()]

          # Calculate area of current piece
          areas[i-1] = gfns.calc_area([x for x, y in indpairs],
                                      [y for x, y in indpairs],
                                                  anl.dx,anl.dy)
            
#==================================================================
# Calculate solidity of pieces (for fragmentation)
#==================================================================

          if anl.addfragmentation or \
             anl.addaxesshape:

            # Define the corners of the pixels
            cornersf = sfns.find_corners(indpairs,
                                      anl.dx,anl.dy)

            # Fit a convex hull to the data points
            hull = ConvexHull(cornersf)
            verts = [(v[0],v[1]) for v in 
                     hull.points[hull.vertices]]

            # Find points in shape
            indpairsCH = sfns.points_in_shape(verts,points)

            # Calculate ratio of object area to shape area
            solid[i-1] = areas[i-1]/gfns.calc_area(
             [x[0] for x in indpairsCH],
             [x[1] for x in indpairsCH],anl.dx,anl.dy)

#==================================================================
# Calculate distance from centroid for all pieces (for dispersion)
#==================================================================

          if (anl.adddispersion or anl.addaxesshape) \
             and numL>1:

            # Define the centroids of the pieces
            centersx[i-1] = np.mean([x for x,y in indpairs])
            centersy[i-1] = np.mean([y for x,y in indpairs])

        # Calculate dispersion
        if (anl.adddispersion or anl.addaxesshape) \
          and numL>1:

          # Calculate distance to center
          dists = np.zeros([numL,numL])
          for i in range(0,numL):
            for j in range(0,numL): 
              dists[i-1,j-1] = gfns.calc_distance(centersx[i-1],
               centersy[i-1],centersx[j-1],centersy[j-1])
          dists = np.ma.array(dists,mask=np.where(dists==0,1,0))
          mdists = np.amin(dists,axis=0)

          areafac = [a/sum(areas) for a in areas]
          distfac = [d/(4*np.sqrt(sum(areas)/math.pi)) \
           for d in mdists]
          dispersion[c] = sum([a*b for a,b in 
                               list(zip(areafac,distfac))])

        elif (anl.adddispersion or anl.addaxesshape) \
          and numL==1: dispersion[c]=0

#==================================================================
# Perimeter and asymmetry prepwork
#==================================================================

        if anl.addperimeter or \
           anl.addasymmetry:

          # Find largest piece and all corner coordinates
          largelabel  = areas.index(max(areas))+1
          indpairslrg = [(round(float(z[1]),2),round(float(z[0]),2))
             for z in df1[df1==largelabel].stack().index.tolist()]
          corners = sfns.find_corners(indpairslrg,anl.dx,anl.dy)

          # Find largest piece and all coordinates
          lonslatslrg = [[x for x, y in corners],
                         [y for x, y in corners]]

#==================================================================
# Add perimeter - only for largest piece
# Uses alpha shapes to define the concave hull with alpha=10 
#==================================================================

        if anl.addperimeter:
        
          # Define the perimeter coordinates
          import alphashape
          hull = alphashape.alphashape(corners, 10)
          hull_pts = hull.exterior.coords.xy

          # Calculate perimeter by summing lengths of all edges
          # Divide by 1000 to get result in km
          perimeter_lp[c] = sum([gfns.calc_distance(
                               hull_pts[0][i],hull_pts[1][i],
                               hull_pts[0][i+1],hull_pts[1][i+1])
                  for i in range(0,len(hull_pts[0])-1)])/1000

#==================================================================
# Add asymmetry
#==================================================================

        if anl.addasymmetry:

          # Fit a convex hull to the data points
          hull = ConvexHull(corners)
          verts = [(v[0],v[1]) for v in 
                     hull.points[hull.vertices]]  

          # Calculate perimeter
          verts.append(verts[0])
          perim = sum([gfns.calc_distance(verts[v][0] ,verts[v][1],
                                       verts[v+1][0],verts[v+1][1]) 
                        for v in range(0,len(verts)-1)])

          # Define Asymmetry
          asymmetry_lp[c] = (4*np.pi*areas[largelabel-1])/\
                                (perim**2)

#==================================================================
# Add fragmentation
#==================================================================

        if anl.addfragmentation or anl.addaxesshape:

          # Calculate metrics related to fragmentation
          connectivity[c]  = 1.-(numL-1.)/(numL+np.log10(
           gfns.calc_area(lonsnzk,latsnzk,anl.dx,anl.dy)))
          solidity[c] = np.mean(solid)
          fragmentation[c] = 1. - (solidity[c]*connectivity[c])
     
#==================================================================
# Add major/minor axes shape
#==================================================================

        if anl.addaxesshape:

          #print("d="+str(dispersion[c]))
          #print("f="+str(fragmentation[c]))
        
          # Only assign shape if f and d are low
          if fragmentation[c]<float(anl.minshapefrag) and \
                dispersion[c]<float(anl.minshapedisp):
            axang[c,:],axlen[c,:] = \
             sfns.fit_ellipse_svd_earth(lonsnzk,latsnzk,
             [dataclon[c],dataclat[c]],fit=True,plot=True)[1:3]

            ellipticity[c] = 1-(axlen[c,1]/axlen[c,0])

#==================================================================
# Convective shape prepwork
#==================================================================

    if anl.addasymmetryc or \
       anl.addfragmentationc or \
       anl.addaxesshapec or \
       anl.adddispersionc or \
       anl.addpiecesc:

      # Find convective locations
      lonsc = lons[k][instrain[k]>anl.convrainthold]
      latsc = lats[k][instrain[k]>anl.convrainthold]
      rainc = instrain[k][instrain[k]>anl.convrainthold]

      if len(latsc)>float(anl.minshapesizec):

        # Define center of convection
        center_c[c,0]= np.mean(lonsc)
        center_c[c,1]= np.mean(latsc)

        # Generate dataframe for object
        dfc = mfns.create_2d_dataframe(lonsc,latsc,
         anl.dx,anl.dy,rainc)

        # Find and label all contiguous areas within object
        labelsc, numLc = sfns.label_wdiags(dfc)
        nyc = [float(i) for i in dfc.index]
        nxc = [float(i) for i in dfc.columns]
        dfc1 = pd.DataFrame(labelsc,[str(i) for i in nyc],
                                    [str(i) for i in nxc])

        # assign number of pieces to pieces
        if anl.addpiecesc:
          pieces_c[c] = numLc

        # Fragmentation prepwork
        if anl.addfragmentationc or \
           anl.addaxesshapec:

          # Predefine solidity array
          solidc = [0.]*numLc
          
          # Define points in the grid
          xgc, ygc = np.meshgrid(nxc,nyc)
          pointsc = [[round(float(x),2) for x in y] 
           for y in np.vstack((xgc.flatten(),ygc.flatten())).T]

        # Loop over all separate pieces
        areasc = [0.]*numLc
        centersxc = [0.]*numLc
        centersyc = [0.]*numLc
        for i in range(1,numLc+1):
        
          # Line 2: Identify pairs of indices for current piece
          # Line 1: Convert all pairs to floats
          indpairsc = [(round(float(z[1]),2),round(float(z[0]),2))
           for z in dfc1[dfc1==i].stack().index.tolist()]

          # Calculate area of current piece
          areasc[i-1] = gfns.calc_area([x for x, y in indpairsc],
                                       [y for x, y in indpairsc],
                                                     anl.dx,anl.dy)
           
#==================================================================
# Calculate solidity of pieces (for fragmentation)
#==================================================================

          if anl.addfragmentationc or \
             anl.addaxesshapec:

            # Define the corners of the pixels
            cornersfc = sfns.find_corners(indpairsc,
                                         anl.dx,anl.dy)

            # Fit a convex hull to the data points
            hullc  = ConvexHull(cornersfc)
            vertsc = [(v[0],v[1]) for v in 
                     hullc.points[hullc.vertices]]

            # Find points in shape
            indpairsCHc = sfns.points_in_shape(vertsc,pointsc)

            # Calculate ratio of object area to shape area
            solidc[i-1] = areasc[i-1]/gfns.calc_area(
             [x[0] for x in indpairsCHc],
             [x[1] for x in indpairsCHc],anl.dx,anl.dy)

#==================================================================
# Calculate distance from centroid for all pieces (for dispersion)
#==================================================================

          if (anl.adddispersionc or \
              anl.addaxesshapec) and numLc>1:

            # Define the centroids of the pieces
            centersxc[i-1] = np.mean([x for x,y in indpairsc])
            centersyc[i-1] = np.mean([y for x,y in indpairsc])

        # Calculate dispersion
        if (anl.adddispersionc or \
            anl.addaxesshapec) and numLc>1:

          # Calculate distance to center
          distsc = np.zeros([numLc,numLc])
          for i in range(0,numLc):
            for j in range(0,numLc): 
              distsc[i-1,j-1] = gfns.calc_distance(centersxc[i-1],
               centersyc[i-1],centersxc[j-1],centersyc[j-1])
          distsc = np.ma.array(distsc,mask=np.where(distsc==0,1,0))
          mdistsc = np.amin(distsc,axis=0)

          areafacc = [a/sum(areasc) for a in areasc]
          distfacc = [d/(4*np.sqrt(sum(areasc)/math.pi)) \
           for d in mdistsc]
          dispersion_c[c] = sum([a*b for a,b in 
                               list(zip(areafacc,distfacc))])

        elif (anl.adddispersionc or \
              anl.addaxesshapec) and numLc==1:
          dispersion_c[c]=0

#==================================================================
# Perimeter and asymmetry prepwork
#==================================================================

        if anl.addasymmetryc:

          # Find largest piece and all corner coordinates
          largelabelc  = areasc.index(max(areasc))+1
          indpairslrgc = [(round(float(z[1]),2),
                           round(float(z[0]),2))
          for z in dfc1[dfc1==largelabelc].stack().index.tolist()]
          cornersc = sfns.find_corners(indpairslrgc,anl.dx,anl.dy)

          # Find largest piece and all coordinates
          lonslatslrgc = [[x for x, y in cornersc],
                          [y for x, y in cornersc]]

          # Fit a convex hull to the data points
          hullc = ConvexHull(cornersc)
          vertsc = [(v[0],v[1]) for v in 
                     hullc.points[hullc.vertices]]

          # Calculate perimeter
          vertsc.append(vertsc[0])
          perimc = sum([gfns.calc_distance(
           vertsc[v][0],vertsc[v][1],
           vertsc[v+1][0],vertsc[v+1][1]) 
           for v in range(0,len(vertsc)-1)])

          # Define Asymmetry
          asymmetry_lp_c[c] = (4*np.pi*areasc[largelabelc-1])/\
                                (perimc**2)

#==================================================================
# Add fragmentation
#==================================================================

        if anl.addfragmentationc or \
           anl.addaxesshapec:

          # Calculate metrics related to fragmentation
          connectivity_c[c]  = 1.-(numLc-1.)/(numLc+np.log10(
           gfns.calc_area(lonsc,latsc,anl.dx,anl.dy)))
          solidity_c[c] = np.mean(solidc)
          fragmentation_c[c] = 1.-(solidity_c[c]*connectivity_c[c])
     
#==================================================================
# Add major/minor axes shape
#==================================================================

        if anl.addaxesshapec:

          #print("dc="+str(dispersion_c[c]))
          #print("fc="+str(fragmentation_c[c]))
        
          if fragmentation_c[c]<float(anl.minshapefragc) and \
                dispersion_c[c]<float(anl.minshapedispc):
            axang_c[c,:],axlen_c[c,:] = \
             sfns.fit_ellipse_svd_earth(lonsc,latsc)[1:3]

            ellipticity_c[c] = 1-(axlen_c[c,1]/axlen_c[c,0])

          ### Plotting for test
          #  if len(latsnzk)>float(anl.minshapesize):
          #    print("Ellipse: yes")
          #    sfns.plot_pf_ellipse(xm,ym,df.values,center_c[c,:],
          #     axang_c[c,:],axlen_c[c,:],fitc)          
          #
          #else:
          #
          #  if len(latsnzk)>float(anl.minshapesize):
          #    print("Ellipse: no")
          #    sfns.plot_pf(xm,ym,df.values)

#==================================================================
# Add land surface information
#==================================================================

    if anl.addlandinfo:
      
      # Import specific libraries
      import shapely.geometry as sgeom
      
      # Get land shape file
      land = gfns.load_land()

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
      TCinfo = gfns.interp_TC(datadtim[c],fTC)

      if TCinfo["TC"]:

        # Check if center is close to TC
        dist_cPF_cTC[c],cPF_in_TC[c],TCname_cPF[c],TCrad_cPF[c] \
          = gfns.calc_if_TC(dataclon[c],dataclat[c],TCinfo,fTC)

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
            = gfns.calc_if_TC(lons[k][il],lats[k][il],TCinfo,fTC)

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
# Write perimeter to file
#==================================================================

  if anl.addperimeter:

    description = "Perimeter of largest piece within the PF. Calculated using alphashapes."
    mfns.write_var("perimeter_lp","Perimeter",description,
     "time",np.float64,"m",fileout,perimeter_lp)

#==================================================================
# Write asymmetry to file
#==================================================================

  if anl.addasymmetry:

    description = "Asymmetry factor for largest piece within PF. 0 = symmetrical (circle). 1 = Highly asymmetrical, non-circular."
    mfns.write_var("asymmetry_lp","Asymmetry factor",description,
     "time",np.float64,"",fileout,asymmetry_lp)

#==================================================================
# Write asymmetry to file
#==================================================================

  if anl.addasymmetryc:

    description = "Asymmetry factor for largest convective piece within PF. 0 = symmetrical (circle). 1 = Highly asymmetrical, non-circular."
    mfns.write_var("asymmetry_lp_c",
     "Asymmetry factor for convective pixels",description,"time",
     np.float64,"",fileout,asymmetry_lp_c)

#==================================================================
# Write fragmentation to file
#==================================================================

  if anl.addfragmentation:

    description = "Fragmentation factor. 0 = One solid piece. 1 = multiple highly fragmented pieces"
    mfns.write_var("fragmentation","Fragmentation factor",
     description,"time",np.float64,"",fileout,fragmentation)

#==================================================================
# Write fragmentation for convective pixels to file
#==================================================================

  if anl.addfragmentationc:

    description = "Fragmentation factor. for convective pixels. 0 = One solid piece. 1 = multiple highly fragmented pieces"
    mfns.write_var("fragmentation_c",
     "Convective fragmentation factor",description,"time",
     np.float64,"",fileout,fragmentation_c)

#==================================================================
# Write dispersion to file
#==================================================================

  if anl.adddispersion:

    description = "Dispersion factor (0-inf). 0 = One piece (no dispersion). 1 = pieces are by a weighted average more than the radius of a circle equivalent to the area of the system away from the center (highly dispersed)."
    mfns.write_var("dispersion","Dispersion factor",
     description,"time",np.float64,"",fileout,dispersion)

#==================================================================
# Write dispersion for convective pixels to file
#==================================================================

  if anl.adddispersionc:

    description = "Dispersion factor (0-inf) for convective pixels. 0 = One piece (no dispersion). 1 = pieces are by a weighted average more than the radius of a circle equivalent to the area of the system away from the center (highly dispersed)."
    mfns.write_var("dispersion_c","Convective dispersion factor",
     description,"time",np.float64,"",fileout,dispersion_c)

#==================================================================
# Write ellipticity to file
#==================================================================

  if anl.addaxesshape:

    description = "Ellipticity factor for PF. Calculated as 1-(major axis length/minor axis length). 1 = highly elliptical. 0 = spherical."
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
# Write ellipticity to file
#==================================================================

  if anl.addaxesshapec:

    description = "Ellipticity factor for convective pixels. Calculated as 1-(major axis length/minor axis length). 1 = highly elliptical. 0 = spherical."
    mfns.write_var("ellipticity_c","Convective ellipticity factor",
      description,"time",np.float64,"",fileout,ellipticity_c)

#==================================================================
# Write length of axes
#==================================================================

    description = "Length of major axis for convective pixels"
    mfns.write_var("mjrax_length_c","Convective major axis length",
     description,"time",np.float64,"m",fileout,axlen_c[:,0])

    description = "Length of minor axis for convective pixels"
    mfns.write_var("mnrax_length_c","Convective minor axis length",
     description,"time",np.float64,"m",fileout,axlen_c[:,1])

#==================================================================
# Write angle major axis makes from north
#==================================================================

    description = "Angle major axis for convective pixels makes with northward vector"
    mfns.write_var("mjrax_angle_c","Convective major axis angle",
     description,"time",np.float64,"degrees",fileout,axang_c[:,0])

    description = "Angle minor axis for convective pixels makes with northward vector"
    mfns.write_var("mnrax_angle_c","Convective minor axis angle",
     description,"time",np.float64,"degrees",fileout,axang_c[:,1])

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
