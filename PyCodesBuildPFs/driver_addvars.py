def driver_addvars(fn):
  """
  Driver to calculate various variables for a single PF.
  fn corresponds to the line number in filenames_av.txt.
  """

#==================================================================
# Begin script
#==================================================================

  # Import libraries
  import PF_functions as PFfunc
  from netCDF4 import Dataset, MFDataset
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
    mjrax_length_lp= [-999]*len(datakeys)
    mnrax_length_lp= [-999]*len(datakeys)
    mjrax_angle_lp = [-999]*len(datakeys)
    mnrax_angle_lp = [-999]*len(datakeys)

  if nl.addasymmetry=="True":
    asymmetry_lp    = [0]*len(datakeys)

  if nl.addfragmentation=="True":
    fragmentation  = [0]*len(datakeys)
    solidity       = [0]*len(datakeys)
    connectivity   = [0]*len(datakeys)

  if nl.addperimeter=="True":
    perimeter_lp     = [0]*len(datakeys)

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
# Assign ERA5 data
#==================================================================

    if nl.addCPE5=="True":
      
      # Import libraries
      import glob
      import datetime as dt
      from scipy.interpolate import interp1d

      # Convert current time to ERA5 time units for comparison
      timestring = str(k)[0:4]+"-"+str(k)[4:6]+"-"+str(k)[6:8]+\
                   " "+str(k)[8:10]+":"+str(k)[10:12]+":00"
      currtime   = PFfunc.time_since(timestring,
                    fd.variables["time"].units)

      # Select which file the time is within
      allfiles = sorted(glob.glob(nl.dataE5dir+nl.fileCPE5id+"*"))
      for fE5 in allfiles:
        ds = Dataset(fE5)
        if ds.variables["time"][0]<=currtime<= \
           ds.variables["time"][-1]:
          times = list(ds.variables["time"][:])
          break
      
      # Select the index(es) of the relevant time(s) 
      if currtime in times:
        timi = times.index(currtime)
      else:
        timi = PFfunc.k_closest(times,currtime,2)   
            
      # Find lat and lon indices
      if dataclon[c]<0: clon = dataclon[c]+360
      latE5 = ds.variables["latitude"][:]
      lonE5 = ds.variables["longitude"][:]
      lati  = np.squeeze([
               PFfunc.k_closest(latE5,dataclat[c]+nl.hda,1), 
               PFfunc.k_closest(latE5,dataclat[c]-nl.hda,1)])
      loni  = np.squeeze([
               PFfunc.k_closest(lonE5,clon-nl.hda,1),
               PFfunc.k_closest(lonE5,clon+nl.hda,1)])

      # Read in subset of data (interpolate in time if necessary)
      if hasattr(timi,"__len__"):
        CAPEnow = np.array(ds.variables["CAPE"][timi[0]:timi[1]+1,
                   lati[0]:lati[1]+1,loni[0]:loni[1]+1])
        intfn = interp1d([times[timi[0]],times[timi[1]]],
                   CAPEnow,axis=0)
        CAPEnow = intfn(currtime)
      else:
        CAPEnow = np.array(ds.variables["CAPE"][timi,
                   lati[0]:lati[1]+1,loni[0]:loni[1]+1])




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
      ar = PFfunc.calc_area(lonsnzk,latsnzk,
        float(nl.dx),float(nl.dy))

      # Convert to units of km**2 
      area[c] = ar/(1000*2)

#==================================================================
# Add area and volumetric rain rate information
#==================================================================

    if nl.addarea=="True" and \
       nl.addvrr=="True":

      # Calculate area and volumetric rain rate
      ar,vrr = PFfunc.calc_area_and_volrainrate(
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
      vrr = PFfunc.calc_area_and_volrainrate(
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
        convarea[c] = PFfunc.calc_area(
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
        convvrr[c] = PFfunc.calc_area_and_volrainrate(
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
        convarea[c],convvrr[c] = PFfunc.calc_area_and_volrainrate(
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
      localsuntime[c] = PFfunc.calc_local_sun_time(
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
          propspd[c],propdir[c] = PFfunc.calc_propagation(
            str(datadtim[c]*100),  dataclon[c],  dataclat[c],
            str(datadtim[c+1]*100),dataclon[c+1],dataclat[c+1])

        # Backward difference for last time
        elif c==len(dataclon)-1:
          propspd[c],propdir[c] = PFfunc.calc_propagation(
            str(datadtim[c-1]*100),dataclon[c-1],dataclat[c-1],
            str(datadtim[c]*100),  dataclon[c],  dataclat[c])

        # Centered difference for all other times
        else:
          propspd[c],propdir[c] = PFfunc.calc_propagation(
            str(datadtim[c-1]*100),dataclon[c-1],dataclat[c-1],
            str(datadtim[c+1]*100),dataclon[c+1],dataclat[c+1])

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
        df = PFfunc.create_2d_dataframe(lonsnzk,latsnzk,
                                        nl.dx,nl.dy,instrainnzk)

        # Find and label all contiguous areas within object
        labels, numL = PFfunc.label_wdiags(df)
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
          areas[i-1] = PFfunc.calc_area([x for x, y in indpairs],
                                        [y for x, y in indpairs],
                                                      nl.dx,nl.dy)

#==================================================================
# Calculate solidity of pieces (for fragmentation)
#==================================================================

          if nl.addfragmentation=="True":

            # Define the corners of the pixels
            cornersf = PFfunc.find_corners(indpairs,
                                      nl.dx,nl.dy)

            # Fit a convex hull to the data points
            hull = ConvexHull(cornersf)
            verts = [(v[0],v[1]) for v in 
                     hull.points[hull.vertices]]

            # Find points in shape
            indpairsCH = PFfunc.points_in_shape(verts,points)

            # Calculate ratio of object area to shape area
            solid[i-1] = areas[i-1]/PFfunc.calc_area(
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
           for z in df1[df1==i].stack().index.tolist()]
        corners = PFfunc.find_corners(indpairslrg,
                                      nl.dx,nl.dy)

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
        perimeter_lp[c] = sum([PFfunc.calc_dist(
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
        perim = sum([PFfunc.calc_dist(verts[v][0]  ,verts[v][1],
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
         PFfunc.calc_area(lonsnzk,latsnzk,nl.dx,nl.dy)))
        solidity[c] = np.mean(solid)
        fragmentation[c] = 1. - (solidity[c]*connectivity[c])
     
#==================================================================
# Add major/minor axes shape
#==================================================================

      if nl.addaxesshape=="True":

        center_lp[c,:],mjrax_length_lp[c],mjrax_angle_lp[c], \
                     mnrax_length_lp[c],mnrax_angle_lp[c] = \
         PFfunc.calc_mjrmnrax(lonslatslrg[0],lonslatslrg[1])

        # Calculate ellipticity  as ratio of major to minor axes
        ellipticity_lp[c] = 1 - (mnrax_length_lp[c]/
                                 mjrax_length_lp[c])

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
      land = PFfunc.load_land()

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
      TCinfo = PFfunc.interp_TC(datadtim[c],fTC)

      if TCinfo["TC"]:

        # Check if center is close to TC
        dist_cPF_cTC[c],cPF_in_TC[c],TCname_cPF[c],TCrad_cPF[c] \
          = PFfunc.calc_if_TC(dataclat[c],dataclon[c],TCinfo,fTC)

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
            = PFfunc.calc_if_TC(lats[k][il],lons[k][il],TCinfo,fTC)

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
# Write normalized time to file
#==================================================================
  
  if nl.addnormtime=="True":

    description = "A fractional time with 0 the first time of the PF and 1 the last time of the PF"
    PFfunc.write_var("normalizedtime","Normalized time",
      description,"time",np.float64,"",fileout,
      normalizedtime,f)

#==================================================================
# Write local time to file
#==================================================================

  if nl.addlocaltime=="True":

    localsuntime = [int(i) for i in localsuntime]
    description = "Calculated as the UTC time plus an offset based on the longitude. The offset is calculated by multiplying the longitude by 24/360. Note: This is not the actual local time. This should typically only be used to calculate times for the diurnal cycle."
    PFfunc.write_var("localsuntime","Local solar time",
      description,"time",int,"",fileout,localsuntime,f)

#==================================================================
# Write max rain rate to file
#==================================================================

  if nl.addmaxrr=="True":

    description = "Maximum rain rate within PF"
    PFfunc.write_var("maxrainrate","Max rain rate",
      description,"time",np.float64,"mm/hr",fileout,
      maxrainrate,f)

#==================================================================
# Write mean rain rate to file
#==================================================================

  if nl.addmeanrr=="True":

    description = "Mean rain rate within PF excluding pixels with zero rain rate"
    PFfunc.write_var("meanrainrate","Mean rain rate",
      description,"time",np.float64,"mm/hr",fileout,
      meanrainrate,f)

#==================================================================
# Write median rain rate to file
#==================================================================

  if nl.addmedianrr=="True":

    description = "Median rain rate within PF excluding pixels with zero rain rate"
    PFfunc.write_var("medianrainrate","Median rain rate",
      description,"time",np.float64,"mm/hr",fileout,
      medianrainrate,f)

#==================================================================
# Write standard deviation rate to file
#==================================================================

  if nl.addstddevrr=="True":

    description = "Standard deviation of rain rate within PF excluding pixels with zero rain rate"
    PFfunc.write_var("stddevrainrate",
      "Standard deviation of rain rate",description,"time",
      np.float64,"mm/hr",fileout,stddevrainrate,f)

#==================================================================
# Write area to file
#==================================================================

  if nl.addarea=="True":

    description = "Area within PF excluding pixels with zero rain rate"
    PFfunc.write_var("area","Area",description,"time",np.float64,
                 "km^2",fileout,area,f)

#==================================================================
# Write volumetric rain rate to file
#==================================================================

  if nl.addvrr=="True":

    description = "Volumetric rain rate within PF excluding pixels with zero rain rate"
    PFfunc.write_var("volrainrate","Volumetric rain rate",
      description,"time",np.float64,"mm hr^-1 km^2",
      fileout,volrainrate,f)

#==================================================================
# Write propagation to file
#==================================================================

  if nl.addpropagation=="True":

    description = "Calculated as the geodesic distance travelled by centroid divided by time taken"
    PFfunc.write_var("propspd","Propagation speed",
      description,"time",np.float64,"m/s",fileout,
      propspd,f)

    description = "Calculated as direction centroid is moving toward from north (clockwise)"
    PFfunc.write_var("propdir","Propagation direction",
      description,"time",np.float64,"degrees",fileout,
      propdir,f)

#==================================================================
# Write convective information to file
#==================================================================

  if nl.addconvrain=="True":

    fileout.conv_rain_threshold = nl.convrainthold
 
    if nl.addconvarea=="True":
      description = "Area of locations with rain rates greater than convective rain rate threshold"
      PFfunc.write_var("conv_area","Convective area",
        description,"time",np.float64,"",fileout,
        convarea,f)

    if nl.addconvvrr=="True":
      description = "Volumetric rain rate of locations with rain rates greater than convective rain rate threshold"
      PFfunc.write_var("conv_vrr",
        "Convective volumetric rain rate",description,
        "time",np.float64,"",fileout,convvrr,f)

    format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."
    description = "Binary indicating if a PF location has a convective rain rate (1 = convective, 0 = not convective)"
    PFfunc.write_group("is_conv_rain",
      "Location has a convective rain rate",description,
      "",format1,fileout,is_conv_rain,f)

#==================================================================
# Write perimeter to file
#==================================================================

  if nl.addperimeter=="True":

    description = "Perimeter of largest piece within the PF. Calculated using alphashapes."
    PFfunc.write_var("perimeter_lp","Perimeter",description,
      "time",np.float64,"m",fileout,perimeter_lp,f)

#==================================================================
# Write asymmetry to file
#==================================================================

  if nl.addasymmetry=="True":

    description = "Asymmetry factor for largest piece within PF. 0 = symmetrical (circle). 1 = Highly asymmetrical, non-circular."
    PFfunc.write_var("asymmetry_lp","Asymmetry factor",description,
      "time",np.float64,"",fileout,asymmetry_lp,f)

#==================================================================
# Write fragmentation to file
#==================================================================

  if nl.addfragmentation=="True":

    description = "Fragmentation factor. 0 = One solid piece. 1 = multiple highly fragmented pieces"
    PFfunc.write_var("fragmentation","Fragmentation factor",
      description,"time",np.float64,"",fileout,fragmentation,f)

#==================================================================
# Write ellipticity to file
#==================================================================

  if nl.addaxesshape=="True":

    description = "Ellipticity factor for largest piece of PF. Calculated as 1-(major axis length/minor axis length). 1 = highly elliptical. 0 = spherical."
    PFfunc.write_var("ellipticity_lp","Ellipticity factor",
      description,"time",np.float64,"",fileout,ellipticity_lp,f)

#==================================================================
# Write coordinates of centroid of largest piece
#==================================================================

    # Write coordinates centroid of largest piece
    lon_center_lp = center_lp[:,0]
    lat_center_lp = center_lp[:,1]

    description = "Longitude of centroid of largest piece"
    PFfunc.write_var("lon_center_lp",
      "Longitude of centroid of largest piece",description,
      "time",np.float64,"degreesE",fileout,lon_center_lp,f)

    description = "Latitude of centroid of largest piece"
    PFfunc.write_var("lat_center_lp",
      "Latitude of centroid of largest piece",description,
      "time",np.float64,"degreesN",fileout,lat_center_lp,f)

#==================================================================
# Write length of axes
#==================================================================

    description = "Length of major axis"
    PFfunc.write_var("mjrax_length_lp","Major axis length",
     description,"time",np.float64,"m",fileout,mjrax_length_lp,f)

    description = "Length of major axis"
    PFfunc.write_var("mnrax_length_lp","Minor axis length",
     description,"time",np.float64,"m",fileout,mnrax_length_lp,f)

#==================================================================
# Write angle major axis makes from north
#==================================================================

    description = "Angle major axis makes with northward vector"
    PFfunc.write_var("mjrax_angle_lp","Major axis angle",
     description,"time",np.float64,"degrees",fileout,
     mjrax_angle_lp,f)

    description = "Length of major axis"
    PFfunc.write_var("mnrax_angle_lp","Minor axis angle",
     description,"time",np.float64,"degrees",fileout,
     mnrax_angle_lp,f)

#==================================================================
# Write boundary information to file
#==================================================================

  if nl.addboundaryinfo=="True":
    
    description = "1 if any part of PF is within a pixel of the tracking domain boundary. Else 0."
    PFfunc.write_var("touchesdombound",
      "PF touches domain boundary",description,"time",
      int,"",fileout,touchesdombound,f)

#==================================================================
# Write land information to file
#==================================================================

  if nl.addlandinfo=="True":
    description = "1 if center of PF is over land. 0 if not."
    PFfunc.write_var("cPF_over_land",
      "Center of PF over land",
      description,"time",int,"",fileout,cPF_over_land,f)

    format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."
    description = "Binary indicating if a PF location is over land (1 = over land, 0 = not over land)"
    PFfunc.write_group("lPF_over_land",
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
      PFfunc.write_var("cPF_in_TC",
        "Proximity of PF center to TC",description,"time",
        int,"",fileout,cPF_in_TC,f)

      # Write dist_cPF_cTC
      description = "Calculated as the geodesic distance of the PF center to TC center"
      PFfunc.write_var("dist_cPF_cTC",
        "Distance of PF center to TC center",description,
        "time",np.float64,"km",fileout,dist_cPF_cTC,f)

      # Write TCradius
      description = "Calculated by finding the largest radius of a TC from IBtracs and doubling."
      PFfunc.write_var("TCrad_cPF","Radius of TC",
        description,"time",np.float64,"km",fileout,
        TCrad_cPF,f)

      # Write list of TCname_cPF as an attribute of PF file
      fileout.TCname_cPF = TCname_cPF
     
      # Create or open a group for lPF_in_TC
      format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."
      description = "Binary indicating if a PF location is within a TC (1 = within TCradius, 0 = not within TC radius)"
      PFfunc.write_group("lPF_in_TC","Location within TC",
        description,"",format1,fileout,lPF_in_TC,f)

      # Create or open a group for dist_lPF_cTC
      description = "Calculated as geodesic distance to TC center from PF location"
      PFfunc.write_group("dist_lPF_cTC",
        "Distance to TC center from PF location",description,
        "km",format1,fileout,dist_lPF_cTC,f)

      # Create or open a group for TCname_lPF
      description = "Name of TC if given location is within TC"
      PFfunc.write_group("TCname_lPF",
        "Name of TC if given location is within TC",
        description,"",format1,fileout,TCname_lPF,f)

      # Create or open a group for TCrad_lPF
      description = "Calculated by finding the largest radius of a TC from IBtracs and doubling."
      PFfunc.write_group("TCrad_lPF","Radius of TC",
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
