def driver_addvars(fn):
  """
  Driver to calculate various variables for a single PF.
  fn corresponds to the line number in filenames_av.txt.
  """

#==================================================================
# Begin script
#==================================================================

  # Import libraries
  import TIPS_functions as fns
  from netCDF4 import Dataset
  import numpy as np
  import pandas as pd
  from scipy.spatial import ConvexHull
  from scipy.spatial import KDTree

  # Read in namelist variables
  nl = Dataset("namelist_av.nc","r")

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
  datadtim = fd.variables["datetime"][:]
  datatim  = fd.variables["time"][:]

  # Get all times 
  datakeys = datalat.__dict__.keys()

#==================================================================
# Preallocate arrays
#==================================================================

  # Preallocate arrays
  if nl.addmaxrr=="True":
    maxrainrate  = [0.]*len(datakeys)

  if nl.addmeanrr=="True":
    meanrainrate = [0.]*len(datakeys)

  if nl.addmedianrr=="True":
    medianrainrate = [0.]*len(datakeys)

  if nl.addstddevrr=="True":
    stddevrainrate = [0.]*len(datakeys)

  if nl.addarea=="True":
    area = [0.]*len(datakeys)

  if nl.addvrr=="True":
    volrainrate = [0.]*len(datakeys)

  if nl.addpropagation=="True":
    propspd = [0.]*len(datakeys)
    propdir = [0.]*len(datakeys)

  if nl.addnormtime=="True":
    # Calculate normalized time variable
    np.seterr(divide='ignore', invalid='ignore')
    normalizedtime = [i for i in range(len(datakeys))]/\
     np.float32(len(datakeys)-1)

  if nl.addlocaltime=="True":
    localsuntime = [""]*len(datakeys)

  if nl.addTCinfo=="True":
    dist_cPF_cTC = [0.]*len(datakeys)
    cPF_in_TC    = [0]*len(datakeys)
    TCname_cPF   = [""]*len(datakeys) 
    TCrad_cPF    = [0.]*len(datakeys)
    lPF_in_TC    = {}
    dist_lPF_cTC = {}
    TCname_lPF   = {}
    TCrad_lPF    = {}
    writeTCdata   = False

  if nl.addlandinfo=="True":
    cPF_over_land = [0]*len(datakeys)
    lPF_over_land = {}

  if nl.addboundaryinfo=="True":
    blatN = float(fd.track_dom_latN)
    blatS = float(fd.track_dom_latS)
    blonE = float(fd.track_dom_lonE)
    blonW = float(fd.track_dom_lonW)
    touchesdombound = [0]*len(datakeys)

  if nl.addconvrain=="True":
    is_conv_rain = {}
    convarea = [0]*len(datakeys)
    convvrr  = [0]*len(datakeys)

  if nl.addaxesshape=="True":
    ellipticity_lp = [-999]*len(datakeys)
    center_lp      = np.zeros((len(datakeys),2))
    center_lp[:,:] = -999
    axlen_lp       = np.zeros((len(datakeys),2))
    axlen_lp[:,:]  = -999
    axang_lp       = np.zeros((len(datakeys),2))
    axang_lp[:,:]  = -999

  if nl.addaxesshapec=="True":
    ellipticity_lp_c = [-999]*len(datakeys)
    center_lp_c      = np.zeros((len(datakeys),2))
    center_lp_c[:,:] = -999
    axlen_lp_c       = np.zeros((len(datakeys),2))
    axlen_lp_c[:,:]  = -999
    axang_lp_c       = np.zeros((len(datakeys),2))
    axang_lp_c[:,:]  = -999

  if nl.addasymmetry=="True":
    asymmetry_lp    = [0]*len(datakeys)

  if nl.addfragmentation=="True":
    fragmentation  = [0]*len(datakeys)
    solidity       = [0]*len(datakeys)
    connectivity   = [0]*len(datakeys)

  if nl.addperimeter=="True":
    perimeter_lp = [0]*len(datakeys)
  
  if nl.addCAPEE5=="True" or nl.addTCWVE5=="True" or \
     nl.addCCTOE5=="True" or nl.addSHRFE5=="True" or \
     nl.addSHRBE5=="True" or nl.addSPHFE5=="True" or \
     nl.addSPHFE5=="True":
    lonsE5 = {}
    latsE5 = {}

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
      ar = fns.calc_area(lonsnzk,latsnzk,
        float(nl.dx),float(nl.dy))

      # Convert to units of km**2 
      area[c] = ar/(1000*2)

#==================================================================
# Add area and volumetric rain rate information
#==================================================================

    if nl.addarea=="True" and \
       nl.addvrr=="True":

      # Calculate area and volumetric rain rate
      ar,vrr = fns.calc_area_and_volrainrate(
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
      vrr = fns.calc_area_and_volrainrate(
        lonsnzk,latsnzk,instrainnzk,
        float(nl.dx),float(nl.dy))[1]

      # Convert to units of mm hr**-1 km**2
      volrainrate[c] = vrr/(1000*2)

#==================================================================
# Add convective rain flag
#==================================================================

    if nl.addconvrain=="True" or \
       nl.addconvarea=="True" or \
       nl.addconvvrr=="True":

      is_conv_rain[k] = np.where(instrain[k]>nl.convrainthold,1,0)

#==================================================================
# Add convective rain area
#==================================================================
  
    if nl.addconvarea=="True" and \
       nl.addconvvrr=="False":
      # Calculate area
      if len(instrain[k][instrain[k]>nl.convrainthold])>0:
        convarea[c] = fns.calc_area(
                    lons[k][instrain[k]>nl.convrainthold],
                    lats[k][instrain[k]>nl.convrainthold],
                    float(nl.dx),float(nl.dy))
      else:
        convarea[c]=0

#==================================================================
# Add convective rain volumetric rain rate
#==================================================================
  
    if nl.addconvarea=="False" and \
       nl.addconvvrr=="True":
      # Calculate area and volumetric rain rate
      if len(instrain[k][instrain[k]>nl.convrainthold])>0:
        convvrr[c] = fns.calc_area_and_volrainrate(
                    lons[k][instrain[k]>nl.convrainthold],
                    lats[k][instrain[k]>nl.convrainthold],
                    instrain[k][instrain[k]>nl.convrainthold],
                    float(nl.dx),float(nl.dy))[1]
      else:
        convvrr[c]=0

#==================================================================
# Add convective rain area and volumetric rain rate
#==================================================================

    if nl.addconvarea=="True" and \
       nl.addconvvrr=="True":
      # Calculate area and volumetric rain rate
      if len(instrain[k][instrain[k]>nl.convrainthold])>0:
        convarea[c],convvrr[c] = fns.calc_area_and_volrainrate(
                    lons[k][instrain[k]>nl.convrainthold],
                    lats[k][instrain[k]>nl.convrainthold],
                    instrain[k][instrain[k]>nl.convrainthold],
                    float(nl.dx),float(nl.dy))
      else:
        convarea[c]=0
        convvrr[c]=0

#==================================================================
# Add local solar time
#==================================================================

    if nl.addlocaltime=="True":
    
      # Calculate local time based on central longitude
      localsuntime[c] = fns.calc_local_sun_time(
                         str(datadtim[c]),dataclon[c])

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

      else:

        # Forward difference for first time
        if c==0:
          propspd[c],propdir[c] = fns.calc_propagation(
            str(datadtim[c]*100),  dataclon[c],  dataclat[c],
            str(datadtim[c+1]*100),dataclon[c+1],dataclat[c+1])

        # Backward difference for last time
        elif c==len(dataclon)-1:
          propspd[c],propdir[c] = fns.calc_propagation(
            str(datadtim[c-1]*100),dataclon[c-1],dataclat[c-1],
            str(datadtim[c]*100),  dataclon[c],  dataclat[c])

        # Centered difference for all other times
        else:
          propspd[c],propdir[c] = fns.calc_propagation(
            str(datadtim[c-1]*100),dataclon[c-1],dataclat[c-1],
            str(datadtim[c+1]*100),dataclon[c+1],dataclat[c+1])

#==================================================================
# Convective shape prepwork
#==================================================================

    # Only calculate if there are enough points and shape is 2d
    if len(latsnzk)>9:

      if nl.addasymmetryc=="True" or \
         nl.addfragmentationc=="True" or \
         nl.addaxesshapec=="True":

        # Generate dataframe for object
        df = fns.create_2d_dataframe(lonsnzk,latsnzk,
                                     nl.dx,nl.dy,instrainnzk)
        #df.loc["2.75"] = 0

        # Find and label all contiguous areas within object
        labels, numL = fns.label_wdiags(df)
        ny = [float(i) for i in df.index]
        nx = [float(i) for i in df.columns]
        df1 = pd.DataFrame(labels,[str(i) for i in ny],
                                  [str(i) for i in nx])

        # Loop over all separate pieces
        areas = [0.]*numL
        for i in range(1,numL+1):
        
          # Line 2: Identify pairs of indices for current piece
          # Line 1: Convert all pairs to floats
          indpairs = [(round(float(z[1]),2),round(float(z[0]),2))
           for z in df1[df1==i].stack().index.tolist()]

          # Calculate area of current piece
          areas[i-1] = fns.calc_area([x for x, y in indpairs],
                        [y for x, y in indpairs],nl.dx,nl.dy)

        # Find convective pixels within largest piece
        largelabel  = areas.index(max(areas))+1
        indpairslrg = [(round(float(z[1]),2),round(float(z[0]),2))
           for z in df1[df1==largelabel].stack().index.tolist()]
        instrainlrg = [df.loc[str(y),str(x)] 
           for x,y in indpairslrg]
        instrainlrgc = [a for a in instrainlrg 
                        if a>=nl.convrainthold]
        indpairslrgc = [b for a,b in zip(instrainlrg,indpairslrg)
                        if a>=nl.convrainthold]
        cornersc = fns.find_corners(indpairslrgc,nl.dx,nl.dy)
        lonslatslrgc = [np.array([x for x, y in cornersc]),
                        np.array([y for x, y in cornersc])]

#==================================================================
# Convective axes shape
#==================================================================

      if nl.addaxesshapec=="True":

        center_lp_c[c,:],axang_lp_c[c,:],axlen_lp_c[c,:] = \
         fns.fit_ellipse_svd(lonslatslrgc[0],lonslatslrgc[1],plot=True)

        ellipticity_lp_c[c] = 1-(axlen_lp_c[c,0]/axlen_lp_c[c,1])


        # Generate dataframe for convective pixels
        #dflrgc = fns.create_2d_dataframe(
        # [x[0] for x in indpairslrgc],[x[1] for x in indpairslrgc],
        # nl.dx,nl.dy,instrainlrgc)

        #import seaborn as sns
        #import matplotlib.pyplot as plt
        #sns.heatmap(df,cmap="GnBu")
        #plt.show()
        #sns.heatmap(dflrgc,cmap="GnBu")
        #plt.show()

#==================================================================
# Shape prepwork
#==================================================================

    # Only calculate if there are enough points and shape is 2d
    if len(latsnzk)>9:

      if nl.addperimeter=="True" or \
         nl.addasymmetry=="True" or \
         nl.addfragmentation=="True" or \
         nl.addaxesshape=="True":

        # Generate dataframe for object
        df = fns.create_2d_dataframe(lonsnzk,latsnzk,
                                     nl.dx,nl.dy,instrainnzk)

        # Find and label all contiguous areas within object
        labels, numL = fns.label_wdiags(df)
        ny = [float(i) for i in df.index]
        nx = [float(i) for i in df.columns]
        df1 = pd.DataFrame(labels,[str(i) for i in ny],
                                  [str(i) for i in nx])

        # Fragmentation prepwork
        if nl.addfragmentation=="True":

          # Predefine solidity array
          solid = [0.]*numL
          
          # Define points in the grid
          xg, yg = np.meshgrid(nx,ny)
          points = [[round(float(x),2) for x in y] 
           for y in np.vstack((xg.flatten(),yg.flatten())).T]

        # Loop over all separate pieces
        areas = [0.]*numL
        for i in range(1,numL+1):
        
          # Line 2: Identify pairs of indices for current piece
          # Line 1: Convert all pairs to floats
          indpairs = [(round(float(z[1]),2),round(float(z[0]),2))
           for z in df1[df1==i].stack().index.tolist()]

          # Calculate area of current piece
          areas[i-1] = fns.calc_area([x for x, y in indpairs],
                                     [y for x, y in indpairs],
                                                  nl.dx,nl.dy)

#==================================================================
# Calculate solidity of pieces (for fragmentation)
#==================================================================

          if nl.addfragmentation=="True":

            # Define the corners of the pixels
            cornersf = fns.find_corners(indpairs,
                                      nl.dx,nl.dy)

            # Fit a convex hull to the data points
            hull = ConvexHull(cornersf)
            verts = [(v[0],v[1]) for v in 
                     hull.points[hull.vertices]]

            # Find points in shape
            indpairsCH = fns.points_in_shape(verts,points)

            # Calculate ratio of object area to shape area
            solid[i-1] = areas[i-1]/fns.calc_area(
             [x[0] for x in indpairsCH],
             [x[1] for x in indpairsCH],nl.dx,nl.dy)

#==================================================================
# Perimeter and asymmetry prepwork
#==================================================================

      if nl.addperimeter=="True" or \
         nl.addasymmetry=="True" or \
         nl.addaxesshape=="True":

        # Find largest piece and all corner coordinates
        largelabel  = areas.index(max(areas))+1
        indpairslrg = [(round(float(z[1]),2),round(float(z[0]),2))
           for z in df1[df1==largelabel].stack().index.tolist()]
        corners = fns.find_corners(indpairslrg,nl.dx,nl.dy)

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
        perimeter_lp[c] = sum([fns.calc_dist(
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
        perim = sum([fns.calc_dist(verts[v][0]  ,verts[v][1],
                                      verts[v+1][0],verts[v+1][1]) 
                          for v in range(0,len(verts)-1)])

        # Define Asymmetry
        asymmetry_lp[c] = (4*np.pi*areas[largelabel-1])/\
                                (perim**2)

#==================================================================
# Add fragmentation
#==================================================================

      if nl.addfragmentation=="True":

        # Calculate metrics related to fragmentation
        connectivity[c]  = 1.-(numL-1.)/(numL+np.log10(
         fns.calc_area(lonsnzk,latsnzk,nl.dx,nl.dy)))
        solidity[c] = np.mean(solid)
        fragmentation[c] = 1. - (solidity[c]*connectivity[c])
     
#==================================================================
# Add major/minor axes shape
#==================================================================

      if nl.addaxesshape=="True":

        center_lp[c,:],axang_lp[c,:],axlen_lp[c,:] = \
         fns.fit_ellipse_svd(lonslatslrg[0],lonslatslrg[1])

        ellipticity_lp_c[c] = 1-(axlen_lp[c,0]/axlen_lp[c,1])

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
      land = fns.load_land()

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
        del(latlonland)   
     
#==================================================================
# Add TC information
#==================================================================

    if nl.addTCinfo=="True":
      
      # Read tropical cyclone data
      fTC = Dataset(nl.dataTCdir+nl.fileTCid) 

      # Interpolate times and latitudes and get radius
      TCinfo = fns.interp_TC(datadtim[c],fTC)

      if TCinfo["TC"]:

        # Check if center is close to TC
        dist_cPF_cTC[c],cPF_in_TC[c],TCname_cPF[c],TCrad_cPF[c] \
          = fns.calc_if_TC(dataclat[c],dataclon[c],TCinfo,fTC)

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
            = fns.calc_if_TC(lats[k][il],lons[k][il],TCinfo,fTC)

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
# ERA5 data prepwork
#==================================================================

    if nl.addCAPEE5=="True" or nl.addTCWVE5=="True" or \
       nl.addCCTOE5=="True" or nl.addSHRFE5=="True" or \
       nl.addSHRBE5=="True" or nl.addSPHFE5=="True" or \
       nl.addSPHFE5=="True":

      import ERA5_functions as E5fns

      # Find file and time indices
      if nl.addCAPEE5=="True":
        fileid1 = nl.fileCAPEE5id
      elif nl.addTCWVE5=="True":
        fileid1 = nl.fileTCWVE5id
      elif nl.addCCTOE5=="True":
        fileid1 = nl.fileCCTOE5id
      elif nl.addSHRFE5=="True":
        fileid1 = nl.fileUSRFE5id
      elif nl.addSHRBE5=="True":
        fileid1 = nl.fileUSRBE5id
      elif nl.addSPHFE5=="True":
        fileid1 = nl.fileSPHFE5id
      elif nl.addSPHBE5=="True":
        fileid1 = nl.fileSPHBE5id
      timestr = str(k)[0:4]+"-"+str(k)[4:6]+"-"+str(k)[6:8]+\
                " "+str(k)[8:10]+":"+str(k)[10:12]+":00"
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,fileid1,timestr)

      # Find coordinates and indices
      loni,lati,lonsE5[k],latsE5[k] = \
       E5fns.get_E5_ss_2D_coords(
        fh,dataclon[c],dataclat[c],nl.hda)
      fh.close()

      xE5 = np.linspace(-nl.hda,nl.hda,len(lonsE5[k]))
      yE5 = np.linspace(-nl.hda,nl.hda,len(latsE5[k]))

#==================================================================
# Assign ERA5 CAPE data
#==================================================================

    if nl.addCAPEE5=="True":

      # Preallocate array
      if c==0:
        CAPEE5 = np.zeros((len(datakeys),len(yE5),len(xE5)))

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileCAPEE5id,timestr)
      
      # Get a subset of the CAPE data
      CAPEE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fh,"CAPE",timi,loni,lati,times,ctime)
      CAPEE5units = fh.variables["CAPE"].units

#==================================================================
# Assign ERA5 TCWV data
#==================================================================

    if nl.addTCWVE5=="True":

      # Preallocate array
      if c==0:
        TCWVE5 = np.zeros((len(datakeys),len(yE5),len(xE5)))

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileTCWVE5id,timestr)
      
      # Get a subset of the CAPE data
      TCWVE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fh,"TCWV",timi,loni,lati,times,ctime)
      TCWVE5units = fh.variables["TCWV"].units

#==================================================================
# Assign ERA5 SPHF data
#==================================================================

    if nl.addSPHFE5=="True":

      # Preallocate array
      if c==0:
        SPHFE5 = np.zeros((len(datakeys),len(yE5),len(xE5)))

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileSPHFE5id,timestr)
      
      # Get a subset of the SPHF data
      SPHFE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fh,"Q",timi,loni,lati,times,ctime)
      SPHFE5units = fh.variables["Q"].units

#==================================================================
# Assign ERA5 SPHB data
#==================================================================

    if nl.addSPHBE5=="True":

      # Preallocate array
      if c==0:
        SPHBE5 = np.zeros((len(datakeys),len(yE5),len(xE5)))

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileSPHBE5id,timestr)
      
      # Get a subset of the SPHF data
      SPHBE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fh,"Q",timi,loni,lati,times,ctime)
      SPHBE5units = fh.variables["Q"].units

#==================================================================
# Assign ERA5 SHRF data
#==================================================================

    if nl.addSHRFE5=="True":

      # Preallocate array
      if c==0:
        USRFE5 = np.zeros((len(datakeys),len(yE5),len(xE5)))
        VSRFE5 = np.zeros((len(datakeys),len(yE5),len(xE5)))
      
      # Find file and time indices
      fhU,timiU,timesU,ctimeU = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileUSRFE5id,timestr)
      fhV,timiV,timesV,ctimeV = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileVSRFE5id,timestr)

      # Get a subset of the SPHF data
      USRFE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fhU,"USHR",timiU,loni,lati,timesU,ctimeU)
      VSRFE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fhV,"VSHR",timiV,loni,lati,timesV,ctimeV)
      SHRFE5units = fhU.variables["USHR"].units
      fhU.close()
      fhV.close()

#==================================================================
# Assign ERA5 SHRB data
#==================================================================

    if nl.addSHRBE5=="True":

      # Preallocate array
      if c==0:
        USRBE5 = np.zeros((len(datakeys),len(yE5),len(xE5)))
        VSRBE5 = np.zeros((len(datakeys),len(yE5),len(xE5)))
      
      # Find file and time indices
      fhU,timiU,timesU,ctimeU = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileUSRBE5id,timestr)
      fhV,timiV,timesV,ctimeV = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileVSRBE5id,timestr)

      # Get a subset of the SPHB data
      USRBE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fhU,"USHR",timiU,loni,lati,timesU,ctimeU)
      VSRBE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fhV,"VSHR",timiV,loni,lati,timesV,ctimeV)
      SHRBE5units = fhU.variables["USHR"].units
      fhU.close()
      fhV.close()

#==================================================================
# Assign ERA5 CCTO data
#==================================================================

    if nl.addCCTOE5=="True":

      # Preallocate array
      if c==0:
        CCTOE5 = np.zeros((len(datakeys),len(yE5),len(xE5)))

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileCCTOE5id,timestr)
      
      # Get a subset of the SPHF data
      CCTOE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fh,"TCC",timi,loni,lati,times,ctime)
      CCTOE5units = fh.variables["TCC"].units

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
# Write normalized time to file
#==================================================================
  
  if nl.addnormtime=="True":

    description = "A fractional time with 0 the first time of the PF and 1 the last time of the PF"
    fns.write_var("normalizedtime","Normalized time",
      description,"time",np.float64,"",fileout,
      normalizedtime,f)

#==================================================================
# Write local time to file
#==================================================================

  if nl.addlocaltime=="True":

    localsuntime = [int(i) for i in localsuntime]
    description = "Calculated as the UTC time plus an offset based on the longitude. The offset is calculated by multiplying the longitude by 24/360. Note: This is not the actual local time. This should typically only be used to calculate times for the diurnal cycle."
    fns.write_var("localsuntime","Local solar time",
      description,"time",int,"",fileout,localsuntime,f)

#==================================================================
# Write max rain rate to file
#==================================================================

  if nl.addmaxrr=="True":

    description = "Maximum rain rate within PF"
    fns.write_var("maxrainrate","Max rain rate",
      description,"time",np.float64,"mm/hr",fileout,
      maxrainrate,f)

#==================================================================
# Write mean rain rate to file
#==================================================================

  if nl.addmeanrr=="True":

    description = "Mean rain rate within PF excluding pixels with zero rain rate"
    fns.write_var("meanrainrate","Mean rain rate",
      description,"time",np.float64,"mm/hr",fileout,
      meanrainrate,f)

#==================================================================
# Write median rain rate to file
#==================================================================

  if nl.addmedianrr=="True":

    description = "Median rain rate within PF excluding pixels with zero rain rate"
    fns.write_var("medianrainrate","Median rain rate",
      description,"time",np.float64,"mm/hr",fileout,
      medianrainrate,f)

#==================================================================
# Write standard deviation rate to file
#==================================================================

  if nl.addstddevrr=="True":

    description = "Standard deviation of rain rate within PF excluding pixels with zero rain rate"
    fns.write_var("stddevrainrate",
      "Standard deviation of rain rate",description,"time",
      np.float64,"mm/hr",fileout,stddevrainrate,f)

#==================================================================
# Write area to file
#==================================================================

  if nl.addarea=="True":

    description = "Area within PF excluding pixels with zero rain rate"
    fns.write_var("area","Area",description,"time",np.float64,
                 "km^2",fileout,area,f)

#==================================================================
# Write volumetric rain rate to file
#==================================================================

  if nl.addvrr=="True":

    description = "Volumetric rain rate within PF excluding pixels with zero rain rate"
    fns.write_var("volrainrate","Volumetric rain rate",
      description,"time",np.float64,"mm hr^-1 km^2",
      fileout,volrainrate,f)

#==================================================================
# Write propagation to file
#==================================================================

  if nl.addpropagation=="True":

    description = "Calculated as the geodesic distance travelled by centroid divided by time taken"
    fns.write_var("propspd","Propagation speed",
      description,"time",np.float64,"m/s",fileout,
      propspd,f)

    description = "Calculated as direction centroid is moving toward from north (clockwise)"
    fns.write_var("propdir","Propagation direction",
      description,"time",np.float64,"degrees",fileout,
      propdir,f)

#==================================================================
# Write convective information to file
#==================================================================

  if nl.addconvrain=="True":

    fileout.conv_rain_threshold = nl.convrainthold
 
    if nl.addconvarea=="True":
      description = "Area of locations with rain rates greater than convective rain rate threshold"
      fns.write_var("conv_area","Convective area",
        description,"time",np.float64,"",fileout,
        convarea,f)

    if nl.addconvvrr=="True":
      description = "Volumetric rain rate of locations with rain rates greater than convective rain rate threshold"
      fns.write_var("conv_vrr",
        "Convective volumetric rain rate",description,
        "time",np.float64,"",fileout,convvrr,f)

    format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."
    description = "Binary indicating if a PF location has a convective rain rate (1 = convective, 0 = not convective)"
    fns.write_group("is_conv_rain",
      "Location has a convective rain rate",description,
      "",format1,fileout,is_conv_rain,f)

#==================================================================
# Write perimeter to file
#==================================================================

  if nl.addperimeter=="True":

    description = "Perimeter of largest piece within the PF. Calculated using alphashapes."
    fns.write_var("perimeter_lp","Perimeter",description,
      "time",np.float64,"m",fileout,perimeter_lp,f)

#==================================================================
# Write asymmetry to file
#==================================================================

  if nl.addasymmetry=="True":

    description = "Asymmetry factor for largest piece within PF. 0 = symmetrical (circle). 1 = Highly asymmetrical, non-circular."
    fns.write_var("asymmetry_lp","Asymmetry factor",description,
      "time",np.float64,"",fileout,asymmetry_lp,f)

#==================================================================
# Write fragmentation to file
#==================================================================

  if nl.addfragmentation=="True":

    description = "Fragmentation factor. 0 = One solid piece. 1 = multiple highly fragmented pieces"
    fns.write_var("fragmentation","Fragmentation factor",
      description,"time",np.float64,"",fileout,fragmentation,f)

#==================================================================
# Write ellipticity to file
#==================================================================

  if nl.addaxesshape=="True":

    description = "Ellipticity factor for largest piece of PF. Calculated as 1-(major axis length/minor axis length). 1 = highly elliptical. 0 = spherical."
    fns.write_var("ellipticity_lp","Ellipticity factor",
      description,"time",np.float64,"",fileout,ellipticity_lp,f)

#==================================================================
# Write coordinates of centroid of largest piece
#==================================================================

    # Write coordinates centroid of largest piece
    lon_center_lp = center_lp[:,0]
    lat_center_lp = center_lp[:,1]

    description = "Longitude of centroid of largest piece"
    fns.write_var("lon_center_lp",
      "Longitude of centroid of largest piece",description,
      "time",np.float64,"degreesE",fileout,lon_center_lp,f)

    description = "Latitude of centroid of largest piece"
    fns.write_var("lat_center_lp",
      "Latitude of centroid of largest piece",description,
      "time",np.float64,"degreesN",fileout,lat_center_lp,f)

#==================================================================
# Write length of axes
#==================================================================

    mjrax_length_lp = axlen_lp[:,0]
    mnrax_length_lp = axlen_lp[:,1]

    description = "Length of major axis"
    fns.write_var("mjrax_length_lp","Major axis length",
     description,"time",np.float64,"m",fileout,mjrax_length_lp,f)

    description = "Length of minor axis"
    fns.write_var("mnrax_length_lp","Minor axis length",
     description,"time",np.float64,"m",fileout,mnrax_length_lp,f)

#==================================================================
# Write angle major axis makes from north
#==================================================================

    mjrax_angle_lp = axang_lp[:,0]
    mnrax_angle_lp = axang_lp[:,1]

    description = "Angle major axis makes with northward vector"
    fns.write_var("mjrax_angle_lp","Major axis angle",
     description,"time",np.float64,"degrees",fileout,
     mjrax_angle_lp,f)

    description = "Angle minor axis makes with northward vector"
    fns.write_var("mnrax_angle_lp","Minor axis angle",
     description,"time",np.float64,"degrees",fileout,
     mnrax_angle_lp,f)

#==================================================================
# Write boundary information to file
#==================================================================

  if nl.addboundaryinfo=="True":
    
    description = "1 if any part of PF is within a pixel of the tracking domain boundary. Else 0."
    fns.write_var("touchesdombound",
     "PF touches domain boundary",description,"time",int,"",
      fileout,touchesdombound,f)

#==================================================================
# Write land information to file
#==================================================================

  if nl.addlandinfo=="True":
    description = "1 if center of PF is over land. 0 if not."
    fns.write_var("cPF_over_land","Center of PF over land",
      description,"time",int,"",fileout,cPF_over_land,f)

    format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."
    description = "Binary indicating if a PF location is over land (1 = over land, 0 = not over land)"
    fns.write_group("lPF_over_land",
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
      fns.write_var("cPF_in_TC",
        "Proximity of PF center to TC",description,"time",
        int,"",fileout,cPF_in_TC,f)

      # Write dist_cPF_cTC
      description = "Calculated as the geodesic distance of the PF center to TC center"
      fns.write_var("dist_cPF_cTC",
        "Distance of PF center to TC center",description,
        "time",np.float64,"km",fileout,dist_cPF_cTC,f)

      # Write TCradius
      description = "Calculated by finding the largest radius of a TC from IBtracs and doubling."
      fns.write_var("TCrad_cPF","Radius of TC",
        description,"time",np.float64,"km",fileout,
        TCrad_cPF,f)

      # Write list of TCname_cPF as an attribute of PF file
      fileout.TCname_cPF = TCname_cPF
     
      # Create or open a group for lPF_in_TC
      format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."
      description = "Binary indicating if a PF location is within a TC (1 = within TCradius, 0 = not within TC radius)"
      fns.write_group("lPF_in_TC","Location within TC",
        description,"",format1,fileout,lPF_in_TC,f)

      # Create or open a group for dist_lPF_cTC
      description = "Calculated as geodesic distance to TC center from PF location"
      fns.write_group("dist_lPF_cTC",
        "Distance to TC center from PF location",description,
        "km",format1,fileout,dist_lPF_cTC,f)

      # Create or open a group for TCname_lPF
      description = "Name of TC if given location is within TC"
      fns.write_group("TCname_lPF",
        "Name of TC if given location is within TC",
        description,"",format1,fileout,TCname_lPF,f)

      # Create or open a group for TCrad_lPF
      description = "Calculated by finding the largest radius of a TC from IBtracs and doubling."
      fns.write_group("TCrad_lPF","Radius of TC",
        description,"",format1,fileout,TCrad_lPF,f)

    else:
    # If not, just set TC attribute as False
      fileout.within_TC = "False"

#==================================================================
# Write coordinate data for environment information to file
#==================================================================

  if nl.addCAPEE5=="True" or nl.addTCWVE5=="True" or \
     nl.addCCTOE5=="True" or nl.addSHRFE5=="True" or \
     nl.addSHRBE5=="True" or nl.addSPHFE5=="True" or \
     nl.addSPHFE5=="True":

    format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."

    description = "longitudes corresponding to ERA5 data"
    fns.write_group("lonsE5","ERA5 longitudes",description,
                    "degreesE",format1,fileout,lonsE5,f)

    description = "latitudes corresponding to ERA5 data"
    fns.write_group("latsE5","ERA5 latitudes",description,
                    "degreesN",format1,fileout,latsE5,f)

    try: xE51 = fileout.createDimension('xE5',len(xE5))
    except: print("xE5 already defined")
    try: yE51 = fileout.createDimension('yE5',len(yE5))
    except: print("yE5 already defined")

    description = "Corresponds to ERA5 data centered on precipitation system centroid. Negative is west. Positive is east."
    fns.write_var("xE5","ERA5 zonal distance from centroid",
     description,("xE5"),np.float64,"degrees",fileout,xE5,f)
    description = "Corresponds to ERA5 data centered on precipitation system centroid. Negative is south. Positive is north."
    fns.write_var("yE5","Meridional distance from centroid",
     description,("yE5"),np.float64,"degrees",fileout,yE5,f)

#==================================================================
# Write CAPE information to file
#==================================================================

  if nl.addCAPEE5=="True":

    description = "Surface-based CAPE from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("CAPE_E5",
     "ERA5 Convective Available Potential Energy",description,
     ("time","yE5","xE5"),np.float64,CAPEE5units,fileout,CAPEE5,f)

#==================================================================
# Write TCWV information to file
#==================================================================

  if nl.addTCWVE5=="True":

    description = "Total column water vapor from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("TCWV_E5","ERA5 Total Column Water Vapor",
     description,("time","yE5","xE5"),np.float64,TCWVE5units,
     fileout,TCWVE5,f)

#==================================================================
# Write SPHF information to file
#==================================================================

  if nl.addSPHFE5=="True":

    description = "Specific humidity between 850-200 hPa from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("SPHU_850-200_E5",
     "ERA5 850-200 hPa mean specific humidity",description,
     ("time","yE5","xE5"),np.float64,SPHFE5units,fileout,SPHFE5,f)

#==================================================================
# Write SPHB information to file
#==================================================================

  if nl.addSPHBE5=="True":

    description = "Specific humidity between 1000-850 hPa from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("SPHU_1000-850_E5",
     "ERA5 1000-850 hPa mean specific humidity",description,
     ("time","yE5","xE5"),np.float64,SPHBE5units,fileout,SPHBE5,f)

#==================================================================
# Write SHRF information to file
#==================================================================

  if nl.addSHRFE5=="True":

    description = "Zonal wind shear between 850-200 hPa from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("USHR_850-200_E5",
     "ERA5 850-200 hPa zonal wind shear",description,
     ("time","yE5","xE5"),np.float64,SHRFE5units,fileout,USRFE5,f)

    description = "Meridional wind shear between 850-200 hPa from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("VSHR_850-200_E5",
     "ERA5 850-200 hPa meridional wind shear",description,
     ("time","yE5","xE5"),np.float64,SHRFE5units,fileout,VSRFE5,f)

#==================================================================
# Write SHRB information to file
#==================================================================

  if nl.addSHRBE5=="True":

    description = "Zonal wind shear between 1000-850 hPa from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("USHR_1000-850_E5",
     "ERA5 1000-850 hPa zonal wind shear",description,
     ("time","yE5","xE5"),np.float64,SHRBE5units,fileout,USRBE5,f)

    description = "Meridional wind shear between 1000-850 hPa from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("VSHR_1000-850_E5",
     "ERA5 1000-850 hPa meridional wind shear",description,
     ("time","yE5","xE5"),np.float64,SHRBE5units,fileout,VSRBE5,f)

#==================================================================
# Write CCTO information to file
#==================================================================

  if nl.addCCTOE5=="True":

    description = "Total cloud cover from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("CCTO_E5","ERA5 Total Cloud Cover",
     description,("time","yE5","xE5"),np.float64,CCTOE5units,
     fileout,CCTOE5,f)

#==================================================================
# Close current file
#==================================================================

  fileout.close()

#==================================================================
# End processing of current PF
#==================================================================
