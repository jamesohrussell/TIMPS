# Import libraries
from netCDF4 import Dataset
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
import sys
import math

def driver_addvars(fn):
  """
  Driver to calculate various variables for a single PF.
  fn corresponds to the line number in filenames_av.txt.
  """

#==================================================================
# Begin script
#==================================================================

  # Read in namelist variables
  nl = Dataset("namelist_av.nc","r")

  # Import custom libraries
  sys.path.insert(0,nl.fnsdir)
  import geophys_functions as gfns
  import time_functions as tfns
  import shape_functions as sfns
  import misc_functions as mfns

  # Read in filename for PF
  for i, row in enumerate(open("filenames_av.txt")):
    if i==fn:
      f = row[:-1]

  print("Working with file "+str(f))

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

  if nl.addmaxrr=="True": maxrainrate  = [-999.]*len(datakeys)

  if nl.addmeanrr=="True": meanrainrate = [-999.]*len(datakeys)

  if nl.addmedianrr=="True": medianrainrate = [-999.]*len(datakeys)

  if nl.addstddevrr=="True": stddevrainrate = [-999.]*len(datakeys)

  if nl.addpieces=="True": pieces = [-999]*len(datakeys)
  if nl.addpiecesc=="True": pieces_c = [-999]*len(datakeys)

  if nl.addarea=="True": area = [-999.]*len(datakeys)

  if nl.addvrr=="True": volrainrate = [-999.]*len(datakeys)

  if nl.addpropagation=="True":
    propspd = [-999.]*len(datakeys)
    propdir = [-999.]*len(datakeys)
    propspdw = [-999.]*len(datakeys)
    propdirw = [-999.]*len(datakeys)

  if nl.addlocaltime=="True": localsolartime = [""]*len(datakeys)

  if nl.addTCinfo=="True":
    dist_cPF_cTC = [-999.]*len(datakeys)
    cPF_in_TC    = [-999]*len(datakeys)
    TCname_cPF   = [""]*len(datakeys) 
    TCrad_cPF    = [-999.]*len(datakeys)
    lPF_in_TC    = {}
    dist_lPF_cTC = {}
    TCname_lPF   = {}
    TCrad_lPF    = {}
    writeTCdata   = False

  if nl.addlandinfo=="True":
    cPF_over_land = [-999]*len(datakeys)
    lPF_over_land = {}

  if nl.addboundaryinfo=="True":
    blatN = float(fd.track_dom_latN)
    blatS = float(fd.track_dom_latS)
    blonE = float(fd.track_dom_lonE)
    blonW = float(fd.track_dom_lonW)
    touchesdombound = [-999]*len(datakeys)

  if nl.addconvrain=="True":
    is_conv_rain = {}
    convarea = [-999.]*len(datakeys)
    convvrr  = [-999.]*len(datakeys)

  if nl.addaxesshape=="True":
    ellipticity = [-999.]*len(datakeys)
    axlen       = np.zeros((len(datakeys),2))
    axlen[:,:]  = -999.
    axang       = np.zeros((len(datakeys),2))
    axang[:,:]  = -999.
  if nl.addaxesshapec=="True":
    ellipticity_c = [-999.]*len(datakeys)
    center_c      = np.zeros((len(datakeys),2))
    center_c[:,:] = -999.
    axlen_c       = np.zeros((len(datakeys),2))
    axlen_c[:,:]  = -999.
    axang_c       = np.zeros((len(datakeys),2))
    axang_c[:,:]  = -999.

  if nl.addasymmetry=="True": 
    asymmetry_lp = [-999.]*len(datakeys)
  if nl.addasymmetryc=="True": 
    asymmetry_lp_c = [-999.]*len(datakeys)

  if nl.addfragmentation=="True" or nl.addaxesshape=="True":
    fragmentation  = [-999.]*len(datakeys)
    solidity       = [-999.]*len(datakeys)
    connectivity   = [-999.]*len(datakeys)
  if nl.addfragmentationc=="True" or nl.addaxesshapec=="True":
    fragmentation_c  = [-999.]*len(datakeys)
    solidity_c       = [-999.]*len(datakeys)
    connectivity_c   = [-999.]*len(datakeys)
 
  if nl.adddispersion=="True" or nl.addaxesshape=="True":
    dispersion = [-999.]*len(datakeys)
  if nl.adddispersionc=="True" or nl.addaxesshapec=="True":
    dispersion_c = [-999.]*len(datakeys)

  if nl.addperimeter=="True": perimeter_lp = [-999.]*len(datakeys)

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

    if nl.addmaxrr=="True":

      # Calculate maximum rain rate in PF
      maxrainrate[c] = np.amax(instrainnzk)

#==================================================================
# Calculate mean rain rate
#==================================================================

    if nl.addmeanrr=="True":

      # Calculate mean rain rates but exclude pixels with 
      #  no rain
      meanrainrate[c]   = np.mean(instrainnzk)

#==================================================================
# Calculate median rain rate
#==================================================================

    if nl.addmedianrr=="True":

      # Calculate median rain rates but exclude pixels with 
      #  no rain
      medianrainrate[c] = np.median(instrainnzk)

#==================================================================
# Calculate standard deviation of rain rate
#==================================================================

    if nl.addstddevrr=="True":
      # Calculate standard deviation of rain rates but 
      #  exclude pixels with no rain
      stddevrainrate[c] = np.std(instrainnzk)

#==================================================================
# Add area 
#==================================================================

    if nl.addarea=="True" and \
       nl.addvrr=="False":

      # Calculate area
      ar = gfns.calc_area(lonsnzk,latsnzk,
        float(nl.dx),float(nl.dy))

      # Convert to units of km**2 
      area[c] = ar/(1000*2)

#==================================================================
# Add area and volumetric rain rate information
#==================================================================

    if nl.addarea=="True" and \
       nl.addvrr=="True":

      # Calculate area and volumetric rain rate
      ar,vrr = gfns.calc_area_and_volrainrate(
        lonsnzk,latsnzk,instrainnzk,float(nl.dx),float(nl.dy))

      # Convert units of km**2 and mm hr**-1 km**2
      area[c] = ar/(1000*2)
      volrainrate[c] = vrr/(1000*2)

#==================================================================
# Add area and volumetric rain rate information
#==================================================================

    if nl.addarea=="False" and \
       nl.addvrr=="True":

      # Calculate volumetric rain rate
      vrr = gfns.calc_area_and_volrainrate(
        lonsnzk,latsnzk,instrainnzk,
        float(nl.dx),float(nl.dy))[1]

      # Convert to units of mm hr**-1 km**2
      volrainrate[c] = vrr/(1000*2)

#==================================================================
# Add local solar time
#==================================================================

    if nl.addlocaltime=="True":
    
      # Calculate local time based on central longitude
      localsolartime[c] = tfns.calc_local_solar_time(
       str(datadtim[c])[0:4],str(datadtim[c])[4:6],
       str(datadtim[c])[6:8],str(datadtim[c])[8:10],
       str(datadtim[c])[10:12],str(0),dataclon[c])

#==================================================================
# Add propagation
#==================================================================

    if nl.addpropagation=="True":
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

    if nl.addconvrain=="True":

      is_conv_rain[k] = np.where(instrain[k]>nl.convrainthold,1,0)

#==================================================================
# Add convective rain area
#==================================================================
  
    if nl.addconvarea=="True" and \
       nl.addconvvrr=="False":
      # Calculate area
      if len(instrain[k][instrain[k]>nl.convrainthold])>0:
        convarea[c] = gfns.calc_area(
                    lons[k][instrain[k]>nl.convrainthold],
                    lats[k][instrain[k]>nl.convrainthold],
                    float(nl.dx),float(nl.dy))
      else: convarea[c]=0

#==================================================================
# Add convective rain volumetric rain rate
#==================================================================
  
    if nl.addconvarea=="False" and \
       nl.addconvvrr=="True":
      # Calculate area and volumetric rain rate
      if len(instrain[k][instrain[k]>nl.convrainthold])>0:
        convvrr[c] = gfns.calc_area_and_volrainrate(
                    lons[k][instrain[k]>nl.convrainthold],
                    lats[k][instrain[k]>nl.convrainthold],
                    instrain[k][instrain[k]>nl.convrainthold],
                    float(nl.dx),float(nl.dy))[1]
      else: convvrr[c]=0

#==================================================================
# Add convective rain area and volumetric rain rate
#==================================================================

    if nl.addconvarea=="True" and \
       nl.addconvvrr=="True":
      # Calculate area and volumetric rain rate
      if len(instrain[k][instrain[k]>nl.convrainthold])>0:
        convarea[c],convvrr[c] = gfns.calc_area_and_volrainrate(
                    lons[k][instrain[k]>nl.convrainthold],
                    lats[k][instrain[k]>nl.convrainthold],
                    instrain[k][instrain[k]>nl.convrainthold],
                    float(nl.dx),float(nl.dy))
      else:
        convarea[c]=0
        convvrr[c]=0

#==================================================================
# Start shape code
#==================================================================

    if nl.addperimeter=="True" or \
       nl.addasymmetry=="True" or \
       nl.addfragmentation=="True" or \
       nl.addaxesshape=="True" or \
       nl.adddispersion=="True" or \
       nl.addpieces=="True":

      # Generate dataframe for object
      df = mfns.create_2d_dataframe(lonsnzk,latsnzk,
                                   nl.dx,nl.dy,instrainnzk)

      # Find and label all contiguous areas within object
      labels, numL = sfns.label_wdiags(df)
      ny = [float(i) for i in df.index]
      nx = [float(i) for i in df.columns]
      df1 = pd.DataFrame(labels,[str(i) for i in ny],
                                [str(i) for i in nx])

      # assign number of pieces to pieces
      if nl.addpieces=="True":
        pieces[c] = numL

      # Only calculate if there are enough points
      if len(latsnzk)>float(nl.minshapesize):

        # Fragmentation prepwork
        if nl.addfragmentation=="True" or \
           nl.addaxesshape=="True":

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
                                                  nl.dx,nl.dy)
            
#==================================================================
# Calculate solidity of pieces (for fragmentation)
#==================================================================

          if nl.addfragmentation=="True" or \
             nl.addaxesshape=="True":

            # Define the corners of the pixels
            cornersf = sfns.find_corners(indpairs,
                                      nl.dx,nl.dy)

            # Fit a convex hull to the data points
            hull = ConvexHull(cornersf)
            verts = [(v[0],v[1]) for v in 
                     hull.points[hull.vertices]]

            # Find points in shape
            indpairsCH = sfns.points_in_shape(verts,points)

            # Calculate ratio of object area to shape area
            solid[i-1] = areas[i-1]/gfns.calc_area(
             [x[0] for x in indpairsCH],
             [x[1] for x in indpairsCH],nl.dx,nl.dy)

#==================================================================
# Calculate distance from centroid for all pieces (for dispersion)
#==================================================================

          if (nl.adddispersion=="True" or nl.addaxesshape=="True") \
             and numL>1:

            # Define the centroids of the pieces
            centersx[i-1] = np.mean([x for x,y in indpairs])
            centersy[i-1] = np.mean([y for x,y in indpairs])

        # Calculate dispersion
        if (nl.adddispersion=="True" or nl.addaxesshape=="True") \
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

        elif (nl.adddispersion=="True" or nl.addaxesshape=="True") \
          and numL==1: dispersion[c]=0

#==================================================================
# Perimeter and asymmetry prepwork
#==================================================================

        if nl.addperimeter=="True" or \
           nl.addasymmetry=="True":

          # Find largest piece and all corner coordinates
          largelabel  = areas.index(max(areas))+1
          indpairslrg = [(round(float(z[1]),2),round(float(z[0]),2))
             for z in df1[df1==largelabel].stack().index.tolist()]
          corners = sfns.find_corners(indpairslrg,nl.dx,nl.dy)

          # Find largest piece and all coordinates
          lonslatslrg = [[x for x, y in corners],
                         [y for x, y in corners]]

#==================================================================
# Add perimeter - only for largest piece
# Uses alpha shapes to define the concave hull with alpha=10 
#==================================================================

        if nl.addperimeter=="True":
        
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

        if nl.addasymmetry=="True":

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

        if nl.addfragmentation=="True" or nl.addaxesshape=="True":

          # Calculate metrics related to fragmentation
          connectivity[c]  = 1.-(numL-1.)/(numL+np.log10(
           gfns.calc_area(lonsnzk,latsnzk,nl.dx,nl.dy)))
          solidity[c] = np.mean(solid)
          fragmentation[c] = 1. - (solidity[c]*connectivity[c])
     
#==================================================================
# Add major/minor axes shape
#==================================================================

        if nl.addaxesshape=="True":

          #print("d="+str(dispersion[c]))
          #print("f="+str(fragmentation[c]))
        
          # Only assign shape if f and d are low
          if fragmentation[c]<float(nl.minshapefrag) and \
                dispersion[c]<float(nl.minshapedisp):
            axang[c,:],axlen[c,:] = \
             sfns.fit_ellipse_svd(lonsnzk,latsnzk)[1:3]

            ellipticity[c] = 1-(axlen[c,1]/axlen[c,0])

          ### Plotting for test
          #
          #  xm,ym = np.meshgrid([x-0.05 for x in nx],
          #                      [y-0.05 for y in ny])
          #  sfns.plot_pf_ellipse(xm,ym,df.values,center,
          #   axang[c,:],axlen[c,:],fit)
          #else:
          #
          #  xm,ym = np.meshgrid([x-0.05 for x in nx],
          #                      [y-0.05 for y in ny])
          #  sfns.plot_pf(xm,ym,df.values)

#==================================================================
# Convective shape prepwork
#==================================================================

    if nl.addasymmetryc=="True" or \
       nl.addfragmentationc=="True" or \
       nl.addaxesshapec=="True" or \
       nl.adddispersionc=="True" or \
       nl.addpiecesc=="True":

      # Find convective locations
      lonsc = lons[k][instrain[k]>nl.convrainthold]
      latsc = lats[k][instrain[k]>nl.convrainthold]
      rainc = instrain[k][instrain[k]>nl.convrainthold]

      if len(latsc)>float(nl.minshapesizec):

        # Define center of convection
        center_c[c,0]= np.mean(lonsc)
        center_c[c,1]= np.mean(latsc)

        # Generate dataframe for object
        dfc = mfns.create_2d_dataframe(lonsc,latsc,
         nl.dx,nl.dy,rainc)

        # Find and label all contiguous areas within object
        labelsc, numLc = sfns.label_wdiags(dfc)
        nyc = [float(i) for i in dfc.index]
        nxc = [float(i) for i in dfc.columns]
        dfc1 = pd.DataFrame(labelsc,[str(i) for i in nyc],
                                    [str(i) for i in nxc])

        # assign number of pieces to pieces
        if nl.addpiecesc=="True":
          pieces_c[c] = numLc

        # Fragmentation prepwork
        if nl.addfragmentationc=="True" or \
           nl.addaxesshapec=="True":

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
                                                     nl.dx,nl.dy)
           
#==================================================================
# Calculate solidity of pieces (for fragmentation)
#==================================================================

          if nl.addfragmentationc=="True" or \
             nl.addaxesshapec=="True":

            # Define the corners of the pixels
            cornersfc = sfns.find_corners(indpairsc,
                                         nl.dx,nl.dy)

            # Fit a convex hull to the data points
            hullc  = ConvexHull(cornersfc)
            vertsc = [(v[0],v[1]) for v in 
                     hullc.points[hullc.vertices]]

            # Find points in shape
            indpairsCHc = sfns.points_in_shape(vertsc,pointsc)

            # Calculate ratio of object area to shape area
            solidc[i-1] = areasc[i-1]/gfns.calc_area(
             [x[0] for x in indpairsCHc],
             [x[1] for x in indpairsCHc],nl.dx,nl.dy)

#==================================================================
# Calculate distance from centroid for all pieces (for dispersion)
#==================================================================

          if (nl.adddispersionc=="True" or \
              nl.addaxesshapec=="True") and numLc>1:

            # Define the centroids of the pieces
            centersxc[i-1] = np.mean([x for x,y in indpairsc])
            centersyc[i-1] = np.mean([y for x,y in indpairsc])

        # Calculate dispersion
        if (nl.adddispersionc=="True" or \
            nl.addaxesshapec=="True") and numLc>1:

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

        elif (nl.adddispersionc=="True" or \
              nl.addaxesshapec=="True") and numLc==1:
          dispersion_c[c]=0

#==================================================================
# Perimeter and asymmetry prepwork
#==================================================================

        if nl.addasymmetryc=="True":

          # Find largest piece and all corner coordinates
          largelabelc  = areasc.index(max(areasc))+1
          indpairslrgc = [(round(float(z[1]),2),
                           round(float(z[0]),2))
          for z in dfc1[dfc1==largelabelc].stack().index.tolist()]
          cornersc = sfns.find_corners(indpairslrgc,nl.dx,nl.dy)

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

        if nl.addfragmentationc=="True" or \
           nl.addaxesshapec=="True":

          # Calculate metrics related to fragmentation
          connectivity_c[c]  = 1.-(numLc-1.)/(numLc+np.log10(
           gfns.calc_area(lonsc,latsc,nl.dx,nl.dy)))
          solidity_c[c] = np.mean(solidc)
          fragmentation_c[c] = 1.-(solidity_c[c]*connectivity_c[c])
     
#==================================================================
# Add major/minor axes shape
#==================================================================

        if nl.addaxesshapec=="True":

          #print("dc="+str(dispersion_c[c]))
          #print("fc="+str(fragmentation_c[c]))
        
          if fragmentation_c[c]<float(nl.minshapefragc) and \
                dispersion_c[c]<float(nl.minshapedispc):
            axang_c[c,:],axlen_c[c,:] = \
             sfns.fit_ellipse_svd(lonsc,latsc)[1:3]

            ellipticity_c[c] = 1-(axlen_c[c,1]/axlen_c[c,0])

          ### Plotting for test
          #  if len(latsnzk)>float(nl.minshapesize):
          #    print("Ellipse: yes")
          #    sfns.plot_pf_ellipse(xm,ym,df.values,center_c[c,:],
          #     axang_c[c,:],axlen_c[c,:],fitc)          
          #
          #else:
          #
          #  if len(latsnzk)>float(nl.minshapesize):
          #    print("Ellipse: no")
          #    sfns.plot_pf(xm,ym,df.values)

#==================================================================
# Add boundary information
#==================================================================

    if nl.addboundaryinfo=="True":
      # Check whether PF is on the tracking domain boundary

      # If not scalar
      if hasattr(lats[k], "__len__"):
        # If any part within a pixel of the northern boundary
        if any(lats[k]>=(blatN-(2*float(nl.dy)))):
          touchesdombound[c] = 1
        # If not, check next boundary and so on for all boundaries
        else:
          if any(lats[k]<=(blatS)):
            touchesdombound[c] = 1
          else:
            if any(lons[k]>=(blonE-(2*float(nl.dx)))):
              touchesdombound[c] = 1
            else:
              if any(lons[k]<=(blonW)):
                touchesdombound[c] = 1
      else:
        if lats[k]>=(blatN-(2*float(nl.dy))):
          touchesdombound[c] = 1
        # If not, check next boundary and so on for all boundaries
        else:
          if lats[k]<=(blatS):
            touchesdombound[c] = 1
          else:
            if lons[k]>=(blonE-(2*float(nl.dx))):
              touchesdombound[c] = 1
            else:
              if lons[k]<=(blonW):
                touchesdombound[c] = 1

#==================================================================
# Add land surface information
#==================================================================

    if nl.addlandinfo=="True":
      
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

    if nl.addTCinfo=="True":
      
      # Read tropical cyclone data
      fTC = Dataset(nl.dataTCdir+nl.fileTCid) 

      # Interpolate times and latitudes and get radius
      TCinfo = gfns.interp_TC(datadtim[c],fTC)

      if TCinfo["TC"]:

        # Check if center is close to TC
        dist_cPF_cTC[c],cPF_in_TC[c],TCname_cPF[c],TCrad_cPF[c] \
          = gfns.calc_if_TC(dataclat[c],dataclon[c],TCinfo,fTC)

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
            = gfns.calc_if_TC(lats[k][il],lons[k][il],TCinfo,fTC)

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

  if nl.addlocaltime=="True":

    localsolartime = [int(i) for i in localsolartime]
    description = "Calculated as the UTC time plus an offset based on the longitude. The offset is calculated by multiplying the longitude by 24/360. Note: This is not the actual local time. This should typically only be used to calculate times for the diurnal cycle."
    mfns.write_var("localsolartime","Local solar time",
      description,"time",np.int64,"",fileout,localsolartime,
      f,int(-999))

#==================================================================
# Write max rain rate to file
#==================================================================

  if nl.addmaxrr=="True":

    description = "Maximum rain rate within PF"
    mfns.write_var("maxrainrate","Max rain rate",description,
     "time",np.float64,"mm/hr",fileout,maxrainrate,f,float(-999))

#==================================================================
# Write mean rain rate to file
#==================================================================

  if nl.addmeanrr=="True":

    description = "Mean rain rate within PF excluding pixels with zero rain rate"
    mfns.write_var("meanrainrate","Mean rain rate",description,
     "time",np.float64,"mm/hr",fileout,meanrainrate,f,float(-999))

#==================================================================
# Write median rain rate to file
#==================================================================

  if nl.addmedianrr=="True":

    description = "Median rain rate within PF excluding pixels with zero rain rate"
    mfns.write_var("medianrainrate","Median rain rate",description,
     "time",np.float64,"mm/hr",fileout,medianrainrate,f,
     float(-999))

#==================================================================
# Write standard deviation rate to file
#==================================================================

  if nl.addstddevrr=="True":

    description = "Standard deviation of rain rate within PF excluding pixels with zero rain rate"
    mfns.write_var("stddevrainrate",
     "Standard deviation of rain rate",description,"time",
      np.float64,"mm/hr",fileout,stddevrainrate,f,float(-999))

#==================================================================
# Write pieces to file
#==================================================================

  if nl.addpieces=="True":

    description = "Number of disconnected pieces making up the precipitation system"
    mfns.write_var("pieces","Pieces",description,"time",np.int64,
     "",fileout,pieces,f,int(-999))

#==================================================================
# Write convective pieces to file
#==================================================================

  if nl.addpiecesc=="True":

    description = "Number of disconnected pieces making up the convective component of the system"
    mfns.write_var("pieces_c","Convective pieces",description,
     "time",np.int64,"",fileout,pieces_c,f,int(-999))

#==================================================================
# Write area to file
#==================================================================

  if nl.addarea=="True":

    description = "Area within PF excluding pixels with zero rain rate"
    mfns.write_var("area","Area",description,"time",np.float64,
     "km^2",fileout,area,f,float(-999))

#==================================================================
# Write volumetric rain rate to file
#==================================================================

  if nl.addvrr=="True":

    description = "Volumetric rain rate within PF excluding pixels with zero rain rate"
    mfns.write_var("volrainrate","Volumetric rain rate",
     description,"time",np.float64,"mm hr^-1 km^2",fileout,
     volrainrate,f,float(-999))

#==================================================================
# Write propagation to file
#==================================================================

  if nl.addpropagation=="True":

    description = "Calculated as the geodesic distance travelled by centroid divided by time"
    mfns.write_var("propspd","Propagation speed",description,
     "time",np.float64,"m/s",fileout,propspd,f,float(-999))

    description = "Calculated as direction centroid is moving toward from north (clockwise)"
    mfns.write_var("propdir","Propagation direction",description,
     "time",np.float64,"degrees",fileout,propdir,f,float(-999))

    description = "Calculated as the geodesic distance travelled by precipitation weighted centroid divided by time"
    mfns.write_var("propspdw","Propagation speed (weighted)",
     description,"time",np.float64,"m/s",fileout,propspdw,f,
     float(-999))

    description = "Calculated as direction precipitation weighted centroid is moving toward from north (clockwise)"
    mfns.write_var("propdirw","Propagation direction (weighted)",
     description,"time",np.float64,"degrees",fileout,propdirw,f,
     float(-999))

#==================================================================
# Write convective information to file
#==================================================================

  if nl.addconvrain=="True" or nl.addconvarea=="True" or \
    nl.addconvvrr=="True":
    fileout.conv_rain_threshold = nl.convrainthold
 
  if nl.addconvarea=="True":
    description = "Area of locations with rain rates greater than convective rain rate threshold"
    mfns.write_var("area_c","Convective area",
     description,"time",np.float64,"",fileout,convarea,f,
     float(-999))

  if nl.addconvvrr=="True":
    description = "Volumetric rain rate of locations with rain rates greater than convective rain rate threshold"
    mfns.write_var("vrr_c","Convective volumetric rain rate",
    description,"time",np.float64,"",fileout,convvrr,f,
    float(-999))

  if nl.addconvrain=="True":
    format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."
    description = "Binary indicating if a PF location has a convective rain rate (1 = convective, 0 = not convective)"
    mfns.write_group("is_conv_rain",
     "Location has a convective rain rate",description,"",
     format1,fileout,is_conv_rain,f)

#==================================================================
# Write perimeter to file
#==================================================================

  if nl.addperimeter=="True":

    description = "Perimeter of largest piece within the PF. Calculated using alphashapes."
    mfns.write_var("perimeter_lp","Perimeter",description,
     "time",np.float64,"m",fileout,perimeter_lp,f,float(-999))

#==================================================================
# Write asymmetry to file
#==================================================================

  if nl.addasymmetry=="True":

    description = "Asymmetry factor for largest piece within PF. 0 = symmetrical (circle). 1 = Highly asymmetrical, non-circular."
    mfns.write_var("asymmetry_lp","Asymmetry factor",description,
     "time",np.float64,"",fileout,asymmetry_lp,f,float(-999))

#==================================================================
# Write asymmetry to file
#==================================================================

  if nl.addasymmetryc=="True":

    description = "Asymmetry factor for largest convective piece within PF. 0 = symmetrical (circle). 1 = Highly asymmetrical, non-circular."
    mfns.write_var("asymmetry_lp_c",
     "Asymmetry factor for convective pixels",description,"time",
     np.float64,"",fileout,asymmetry_lp_c,f,float(-999))

#==================================================================
# Write fragmentation to file
#==================================================================

  if nl.addfragmentation=="True":

    description = "Fragmentation factor. 0 = One solid piece. 1 = multiple highly fragmented pieces"
    mfns.write_var("fragmentation","Fragmentation factor",
     description,"time",np.float64,"",fileout,fragmentation,f,
     float(-999))

#==================================================================
# Write fragmentation for convective pixels to file
#==================================================================

  if nl.addfragmentationc=="True":

    description = "Fragmentation factor. for convective pixels. 0 = One solid piece. 1 = multiple highly fragmented pieces"
    mfns.write_var("fragmentation_c",
     "Convective fragmentation factor",description,"time",
     np.float64,"",fileout,fragmentation_c,f,float(-999))

#==================================================================
# Write dispersion to file
#==================================================================

  if nl.adddispersion=="True":

    description = "Dispersion factor (0-inf). 0 = One piece (no dispersion). 1 = pieces are by a weighted average more than the radius of a circle equivalent to the area of the system away from the center (highly dispersed)."
    mfns.write_var("dispersion","Dispersion factor",
     description,"time",np.float64,"",fileout,dispersion,
     f,float(-999))

#==================================================================
# Write dispersion for convective pixels to file
#==================================================================

  if nl.adddispersionc=="True":

    description = "Dispersion factor (0-inf) for convective pixels. 0 = One piece (no dispersion). 1 = pieces are by a weighted average more than the radius of a circle equivalent to the area of the system away from the center (highly dispersed)."
    mfns.write_var("dispersion_c","Convective dispersion factor",
     description,"time",np.float64,"",fileout,dispersion_c,f,
     float(-999))

#==================================================================
# Write ellipticity to file
#==================================================================

  if nl.addaxesshape=="True":

    description = "Ellipticity factor for PF. Calculated as 1-(major axis length/minor axis length). 1 = highly elliptical. 0 = spherical."
    mfns.write_var("ellipticity","Ellipticity factor",
      description,"time",np.float64,"",fileout,ellipticity,
      f,float(-999))

#==================================================================
# Write length of axes
#==================================================================

    description = "Length of major axis"
    mfns.write_var("mjrax_length","Major axis length",
     description,"time",np.float64,"m",fileout,axlen[:,0],
     f,float(-999))

    description = "Length of minor axis"
    mfns.write_var("mnrax_length","Minor axis length",
     description,"time",np.float64,"m",fileout,axlen[:,1],
     f,float(-999))

#==================================================================
# Write angle major axis makes from north
#==================================================================

    description = "Angle major axis makes with northward vector"
    mfns.write_var("mjrax_angle","Major axis angle",
     description,"time",np.float64,"degrees",fileout,
     axang[:,0],f,float(-999))

    description = "Angle minor axis makes with northward vector"
    mfns.write_var("mnrax_angle","Minor axis angle",
     description,"time",np.float64,"degrees",fileout,
     axang[:,1],f,float(-999))

#==================================================================
# Write ellipticity to file
#==================================================================

  if nl.addaxesshapec=="True":

    description = "Ellipticity factor for convective pixels. Calculated as 1-(major axis length/minor axis length). 1 = highly elliptical. 0 = spherical."
    mfns.write_var("ellipticity_c","Convective ellipticity factor",
      description,"time",np.float64,"",fileout,ellipticity_c,
      f,float(-999))

#==================================================================
# Write length of axes
#==================================================================

    description = "Length of major axis for convective pixels"
    mfns.write_var("mjrax_length_c","Convective major axis length",
     description,"time",np.float64,"m",fileout,axlen_c[:,0],
     f,float(-999))

    description = "Length of minor axis for convective pixels"
    mfns.write_var("mnrax_length_c","Convective minor axis length",
     description,"time",np.float64,"m",fileout,axlen_c[:,1],
     f,float(-999))

#==================================================================
# Write angle major axis makes from north
#==================================================================

    description = "Angle major axis for convective pixels makes with northward vector"
    mfns.write_var("mjrax_angle_c","Convective major axis angle",
     description,"time",np.float64,"degrees",fileout,
     axang_c[:,0],f,float(-999))

    description = "Angle minor axis for convective pixels makes with northward vector"
    mfns.write_var("mnrax_angle_c","Convective minor axis angle",
     description,"time",np.float64,"degrees",fileout,
     axang_c[:,1],f,float(-999))

#==================================================================
# Write boundary information to file
#==================================================================

  if nl.addboundaryinfo=="True":
    
    description = "1 if any part of PF is within a pixel of the tracking domain boundary. Else 0."
    mfns.write_var("touchesdombound",
     "PF touches domain boundary",description,"time",np.int64,"",
      fileout,touchesdombound,f,int(-999))

#==================================================================
# Write land information to file
#==================================================================

  if nl.addlandinfo=="True":
    description = "1 if center of PF is over land. 0 if not."
    mfns.write_var("cPF_over_land","Center of PF over land",
      description,"time",np.int64,"",fileout,cPF_over_land,f,int(-999))

    format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."
    description = "Binary indicating if a PF location is over land (1 = over land, 0 = not over land)"
    mfns.write_group("lPF_over_land",
      "Location over land",description,"",
      format1,fileout,lPF_over_land,f)

#==================================================================
# Write TC information to file
#==================================================================

  if nl.addTCinfo=="True":
    
    if writeTCdata: 

    # If any part of PF IS within TC radius set TC attribute as True and write
    #  other variables describing which parts of PF are within TC radius 
      fileout.within_TC = "True"

      # Write cPF_in_TC
      description = "If PF center is within a maximum radius of TC center this is 1, if not, this is 0."
      mfns.write_var("cPF_in_TC",
        "Proximity of PF center to TC",description,"time",
        np.int64,"",fileout,cPF_in_TC,f,int(-999))

      # Write dist_cPF_cTC
      description = "Calculated as the geodesic distance of the PF center to TC center"
      mfns.write_var("dist_cPF_cTC",
        "Distance of PF center to TC center",description,
        "time",np.float64,"km",fileout,dist_cPF_cTC,f,float(-999))

      # Write TCradius
      description = "Calculated by finding the largest radius of a TC from IBtracs and doubling."
      mfns.write_var("TCrad_cPF","Radius of TC",
        description,"time",np.float64,"km",fileout,
        TCrad_cPF,f,float(-999))

      # Write list of TCname_cPF as an attribute of PF file
      fileout.TCname_cPF = TCname_cPF
     
      # Create or open a group for lPF_in_TC
      format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."
      description = "Binary indicating if a PF location is within a TC (1 = within TCradius, 0 = not within TC radius)"
      mfns.write_group("lPF_in_TC","Location within TC",
        description,"",format1,fileout,lPF_in_TC,f)

      # Create or open a group for dist_lPF_cTC
      description = "Calculated as geodesic distance to TC center from PF location"
      mfns.write_group("dist_lPF_cTC",
        "Distance to TC center from PF location",description,
        "km",format1,fileout,dist_lPF_cTC,f)

      # Create or open a group for TCname_lPF
      description = "Name of TC if given location is within TC"
      mfns.write_group("TCname_lPF",
        "Name of TC if given location is within TC",
        description,"",format1,fileout,TCname_lPF,f)

      # Create or open a group for TCrad_lPF
      description = "Calculated by finding the largest radius of a TC from IBtracs and doubling."
      mfns.write_group("TCrad_lPF","Radius of TC",
        description,"",format1,fileout,TCrad_lPF,f)

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
