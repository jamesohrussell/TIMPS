#==================================================================
# Functions used to add variables to IMERG PFs
# Author: James Russell 2019
#
# Tested with:
# * numpy 1.16.3
# * scipy 1.2.1
# * matplotlib 3.0.3
# * pandas 0.24.2
# * cartopy 0.17.0
# * shapely 1.6.4
# * fiona 1.8.6
# * area 1.1.1
# * pyproj 1.9.6
# * geopy 1.20.0
#
# Description of functions:
#
# * calc_area
#   - Calculates area on earth from a list of pixels 
#      locations
#
# * calc_area_and_volrainrate
#   - Calculates these variables on earth from a list of 
#      pixels locations and their rain rates
#
# * calc_distance
#   - Calculates distance between two points on earth 
#
# * calc_distance_and_angle
#   - Calculates distance between two points on earth 
#      and the angle the vector those points make from 
#      north
#
# * calc_propagation
#   - Calculate propagation speed and direction between 
#      two points on earth
#
# * interp_TC
#   - Reads IBtracs database and interpolates position of 
#      all TCs within 3 hours of input time, to the input 
#      time. Also outputs a list of largest possible TC 
#      radii.
#
# * calc_if_TC
#   - Takes output from interp_TC and a latitude, longitude 
#      location and calculates whether a TC is close to 
#      that location. Outputs information on the TC that is 
#      closest if it is within largest possible TC radius.
#
# * load_land
#   - Loads a land shape file
# 
# * is_land
#   - Takes a latitude and longitude location, reads in a 
#      land area shape file, it checks location against 
#      shape file to return True if location is over land.
#
# * calc_local_sun_time
#   - Takes a date and time, and a longitude, and adds an 
#      offset factor to give a local solar time (i.e. for a
#      diurnal cycle). Not the actual local time.
#
# * create_2d_dataframe
#   - Generates a dataframe from lists of coordinates and 
#      the corresponding data.
#
# * label_wdiags
#   - Replicates the scipy.ndimage.label() function but 
#      connects features at diagonals
#
# * find_corners
#   - Finds all x,y coordinates of all corners given a set of 
#      pixel center coordinates
#
# * points_in_shape
#   - Finds all x,y coordinates within a shape as defined by its
#      vertices
#
# * calc_mjrmnrax
#   - Finds the major and minor axis of a set of points on earth
#
# * write_var
#   - Defines/opens a variable and writes the data to a 
#      netcdf file.
#
# * write_group
#   - Defines/opens a group and writes data to a data group
#      within that netcdf group as attribute and value 
#      pairs. This is a means for writing a python 
#      dictionary to a netcdf file.
#
#==================================================================

#==================================================================
# Calculate area 
#==================================================================

def calc_area(lons,lats,dx,dy):
  """
  Calculates the area on earth given a set of pixels.
  
  Inputs:
  1,2) Lists of latitude and longitude coordinates of pixel 
   centers
  3,4) Scalars of the pixel width (dx) and height (dy). 
  All inputs should be in units of degrees longitude or latitude.
  
  Output is an area in units of m**2.

  Requires area 1.1.1 (https://pypi.org/project/area/).
   Note: At time of writing this was not available 
   through conda and had to be installed manually.
  """

  # Import libraries
  from area import area
  import numpy

  # Loop over all pixels
  for l in range(0,len(lats)):
    # Create geojson format variable for each pixel. 
    geom = {'type': 'Polygon',
            'coordinates': [[[lons[l]-(dx/2.), lats[l]-(dy/2.)], 
                             [lons[l]+(dx/2.), lats[l]-(dy/2.)],
                             [lons[l]+(dx/2.), lats[l]+(dy/2.)],
                             [lons[l]-(dx/2.), lats[l]+(dy/2.)],
                             [lons[l]-(dx/2.), lats[l]-(dy/2.)]]]}

    # Calculate area for each individual pixel
    proj_area = area(geom)

    # Sum area of all pixels
    if l==0:
      area1 = proj_area
    else:
      area1 = area1 + proj_area

  # Return area
  return(area1)



#==================================================================
# Calculate area and volumetric rain rate
#==================================================================

def calc_area_and_volrainrate(lons,lats,rain,dx,dy):
  """
  Calculates the area on earth and the volumetric rain 
   rate of that area given a set of pixels.
  
  Inputs:
   1,2) Lists of latitude and longitude coordinates of pixel 
    centers.
   3) List of corresponding rainfalls for each pixel
   4,5) Scalars of the pixel width (dx) and height (dy). 
   All inputs should be in units of degrees longitude or 
    latitude, with rainfall in mm/hr.
  
  Output is an area in units of m**2 and a volumetric rain 
   rate in units of mm/hr m**2.

  Requires area 1.1.1 (https://pypi.org/project/area/).
   Note: At time of writing this was not available 
   through conda and had to be installed manually.
  """

  # Import libraries
  from area import area

  # Loop over all pixels
  for l in range(0,len(lats)):
    # Create geojson format variable for each pixel. 
    geom = {'type': 'Polygon',
            'coordinates': [[[lons[l]-(dx/2.), lats[l]-(dy/2.)], 
                             [lons[l]+(dx/2.), lats[l]-(dy/2.)],
                             [lons[l]+(dx/2.), lats[l]+(dy/2.)],
                             [lons[l]-(dx/2.), lats[l]+(dy/2.)],
                             [lons[l]-(dx/2.), lats[l]-(dy/2.)]]]}

    # Calculate area for each individual pixel
    proj_area = area(geom)

    # Calculate volumetric rain rate per pixel
    volrainratel = proj_area * rain[l]

    # Add area and volumetric rain rates of all pixels
    if l==0:
      area1 = proj_area
      volrainrate = volrainratel
    else:
      area1 = area1 + proj_area
      volrainrate = volrainrate + volrainratel

  # Return area and volumetric rain rate
  return(area1,volrainrate)



#==================================================================
# Calculate distance between two points on earth
#==================================================================

def calc_dist(lon1,lat1,lon2,lat2):
  """
  Calculates distance on earth between two sets of 
   latitude, longitude coordinates.

  Input:
   1,2) lon1, lat1
   3,4) lon2, lat2

  Output is the distance in units of m.

  Requires pyproj 1.9.6 (conda install -c conda-forge pyproj; 
   https://pypi.org/project/pyproj/) and geopy 1.20.0
   (conda install -c conda-forge geopy; 
   https://pypi.org/project/geopy/).
  """

  # Import libraries
  from geopy.distance import geodesic
  import pyproj

  # Get projection
  geod = pyproj.Geod(ellps='WGS84')

  # Calculate propagation speed
  dist = geodesic((lat1,lon1),(lat2,lon2)).m

  return(dist)



#==================================================================
# Calculate distance and angle of a vector on earth
#==================================================================

def calc_distandangle(lon1,lat1,lon2,lat2):
  """
  Calculates distance on earth between two sets of 
   latitude, longitude coordinates, and then calculates 
   angle that vector makes from north.

  Inputs:
  1,2) lat1, lon1
  3,4) lat2, lon2

  Output is the distance in units of m and angle 
   the vector makes from north.

  Requires pyproj 1.9.6 (conda install -c conda-forge pyproj; 
   https://pypi.org/project/pyproj/) and geopy 1.20.0
   (conda install -c conda-forge geopy; 
   https://pypi.org/project/geopy/).
  """

  # Import libraries
  from geopy.distance import geodesic
  import pyproj

  # Get projection
  geod = pyproj.Geod(ellps='WGS84')

  # Calculate propagation speed
  dist = geodesic((lat1,lon1),(lat2,lon2)).m

  # Calculate propagation direction
  az1, az2, di = geod.inv(str(lon1),str(lat1),
                          str(lon2),str(lat2))

  # Adjust for degrees from north
  if az1<0:
    angle = 360+az1
  else:
    angle = az1

  return(dist,angle)



#==================================================================
# Calculate propagation
#==================================================================

def calc_propagation(date1,lon1,lat1,date2,lon2,lat2):
  """
  Calculates propagation speed and direction on earth 
   given two dates and times with corresponding latitudes 
   and longitudes.

  Inputs:
   1,2,3) date and time 1, lat1, lon1
   4,5,6) date and time 2, lat2, lon2. 
   Dates and times should be in string format YYYYMMDDhhmmss. 
   All others should be floats or integers.

  Output is the speed in units of m/s and angle the 
   propagation vector makes from north.

  Requires pyproj 1.9.6 (conda install -c conda-forge pyproj; 
   https://pypi.org/project/pyproj/) and geopy 1.20.0
   (conda install -c conda-forge geopy; 
   https://pypi.org/project/geopy/).
  """

  # Import libraries
  from geopy.distance import geodesic
  import pyproj
  import datetime as dt

  # Get projection
  geod = pyproj.Geod(ellps='WGS84')

  # Make time objects
  dtob1 = dt.datetime(int(date1[0:4]),int(date1[4:6]),
                      int(date1[6:8]),int(date1[8:10]),
                      int(date1[10:12]),int(date1[12:14]))
  dtob2 = dt.datetime(int(date2[0:4]),int(date2[4:6]),
                      int(date2[6:8]),int(date2[8:10]),
                      int(date2[10:12]),int(date2[12:14]))

  # Calculate propagation speed
  propspd = (geodesic((lat1,lon1),(lat2,lon2)).m)/\
             ((dtob2-dtob1).seconds)

  # Calculate propagation direction
  az1, az2, di = geod.inv(str(lon1),str(lat1),
                          str(lon2),str(lat2))

  # Adjust for degrees from north
  if az1<0:
    propdir = 360+az1
  else:
    propdir = az1

  return(propspd,propdir)



#==================================================================
# Interpolate TC information to a time
#==================================================================

def interp_TC(dtim,fTC):
  """
  Interpolates information on a TC to a given a time. 

  Inputs:
  1) A date and time in YYYYMMDDhhmm string format to interpolate
  2) File handle (fTC) for IBTrACS data.

  Output is  a dictionary of information about the TC after
   interpolating that information to the same time as that 
   input. In that dictionary, it returns lists of TC radii,
   lat and lon position of the TC center, and the 
   corresponding indices of the data in the ibtracs data. If
   TC is not within 1.5 hours of the given time, returns 
   dictionary with one key and value pair TC:False. If within 
   1.5 hours this TC key is True and all other information is 
   provided."

  Requires geopy 1.20.0 (conda install -c conda-forge geopy; 
   https://pypi.org/project/geopy/).
  """

  # Import libraries 
  import datetime
  import time
  from geopy.distance import geodesic 
  import numpy as np

  # Convert current time to TC time units
  #  *Hard coded for IBTrACS default unit of days
  timeunits = fTC.variables["time"].units
  d1   = datetime.datetime.strptime(timeunits[11:21], 
                                    "%Y-%m-%d")
  d2   = datetime.datetime.strptime(str(dtim),
                                    "%Y%m%d%H%M")
  d1ts = time.mktime(d1.timetuple())
  d2ts = time.mktime(d2.timetuple())
  TCtimenow = (d2ts-d1ts)/(3600.*24)

  # Check for TC times within just over 1.5 hours
  #  i.e. since TCs observations are typically 3 hours apart
  TCtimelow = TCtimenow-0.08
  TCtimehig = TCtimenow+0.08
  indTC = np.asarray(np.where(
    (fTC.variables["time"][:]>TCtimelow) 
    & (fTC.variables["time"][:]<TCtimehig))).T.tolist()

  if not indTC:
    TCinfo = {"TC":False}
    return(TCinfo)

  TCinfo = {"TC":True}

  # Preallocate some output arrays
  TCradmax = [0.]*len(indTC)
  TClatnow = [0.]*len(indTC)
  TClonnow = [0.]*len(indTC)

  # Loop over all TCs with a close time
  for iTC in range(len(indTC)):
    # Find TC times and locations and interpolate them to
    #  the time of the PF

    # Instance where time of PF and TC ob are same
    if TCtimenow==fTC.variables["time"][indTC[iTC][0],
                                        indTC[iTC][1]]:
      TClatnow[iTC] = fTC.variables["lat"][indTC[iTC][0],
                                           indTC[iTC][1]]
      TClonnow[iTC] = fTC.variables["lon"][indTC[iTC][0],
                                           indTC[iTC][1]]

    # Instance where time of PF is earlier than first TC ob 
    elif ((TCtimenow<fTC.variables["time"][indTC[iTC][0],
              indTC[iTC][1]]) & (indTC[iTC][1]==0)):
      TClatnow[iTC] = fTC.variables["lat"][indTC[iTC][0],
                                           indTC[iTC][1]]
      TClonnow[iTC] = fTC.variables["lon"][indTC[iTC][0],
                                           indTC[iTC][1]]

    # Instance where time of PF is earlier than TC ob 
    elif TCtimenow<fTC.variables["time"][indTC[iTC][0],
                                         indTC[iTC][1]]:
      TClatnow[iTC] = np.interp(TCtimenow,
        fTC.variables["time"][indTC[iTC][0],
        indTC[iTC][1]-1:indTC[iTC][1]+1],
        fTC.variables["lat"][indTC[iTC][0],
        indTC[iTC][1]-1:indTC[iTC][1]+1])
      TClonnow[iTC] = np.interp(TCtimenow,
        fTC.variables["time"][indTC[iTC][0],
        indTC[iTC][1]-1:indTC[iTC][1]+1],
        fTC.variables["lon"][indTC[iTC][0],
        indTC[iTC][1]-1:indTC[iTC][1]+1])

    # Instance where time of PF is later than first TC ob 
    elif TCtimenow>fTC.variables["time"][indTC[iTC][0],
                                         indTC[iTC][1]]:
      TCtim = fTC.variables["time"][indTC[iTC][0],
        indTC[iTC][1]:indTC[iTC][1]+2]
      TCtim = TCtim[~TCtim.mask]
      if len(TCtim)==2:
        TClat = fTC.variables["lat"][indTC[iTC][0],
          indTC[iTC][1]:indTC[iTC][1]+2]
        TClon = fTC.variables["lon"][indTC[iTC][0],
          indTC[iTC][1]:indTC[iTC][1]+2]
        TClatnow[iTC] = np.interp(TCtimenow,TCtim,TClat)
        TClonnow[iTC] = np.interp(TCtimenow,TCtim,TClon)
      else:
        TClatnow[iTC] = fTC.variables["lat"][indTC[iTC][0],
                                             indTC[iTC][1]]
        TClonnow[iTC] = fTC.variables["lon"][indTC[iTC][0],
                                             indTC[iTC][1]]

    # Read in the widest TC radius possible
    TCrad = np.ma.empty(6)
    TCrad[0] = fTC.variables["td9635_roci"][indTC[iTC][0],
                                            indTC[iTC][1]]
    TCrad[1] = fTC.variables["bom_roci"][indTC[iTC][0],
                                         indTC[iTC][1]]
    TCrad[2] = max(fTC.variables["bom_r34"][indTC[iTC][0],
                                            indTC[iTC][1],:])
    TCrad[3] = max(fTC.variables["reunion_r34"][indTC[iTC][0],
                                                indTC[iTC][1],:])
    TCrad[4] = fTC.variables["tokyo_r30_long"][indTC[iTC][0],
                                               indTC[iTC][1]]
    TCrad[5] = max(fTC.variables["usa_r34"][indTC[iTC][0],
                                            indTC[iTC][1],:])

    if TCrad.count()==0:
    # If there is no information on a minimum radius select a 
    #  typical average radius of 0 m/s winds(source: Chavas et
    #  al. 2016: Observed Tropical Cyclone Size Revisited. JCLI)
      TCradmax[iTC] = 600.
    else:
    # If there is information, select the largest and 
    #  convert for nautical miles to km and double
    #  to ensure all of TC is accounted for (i.e. lighter 
    #  than 30 kt winds. Source above shows this is likely 
    #  underestimate)
      TCradmax[iTC] = (np.ma.MaskedArray.max(TCrad)*1.852)*2

  TCinfo["TCrad"] = TCradmax
  TCinfo["TClat"] = TClatnow
  TCinfo["TClon"] = TClonnow
  TCinfo["indTC"] = indTC

  return(TCinfo)



#==================================================================
# Calculate TC information
#==================================================================

def calc_if_TC(lon,lat,TCinfo,fTC):
  """
  Calculates proximity to TC. 

  Inputs:
   1,2) lat and lon locations for position to check
   3) A dictionary with:
    TClat and TClon - position of TC center
    TCrad - list of radiuses 
   4) Open file handle of the IBTrACS dataset

  Outputs are the distance between the lat, lon location 
   and the active TC centers at that time, a boolean 
   indicating if it is within a radius of the TC, the name 
   of the TC it is within radius of, and the radius used to
   estimate if location is within the TC.

  Requires geopy 1.20.0 (conda install -c conda-forge geopy; 
   https://pypi.org/project/geopy/).
  """

  # Import libraries
  from geopy.distance import geodesic 
  import numpy as np
  import netCDF4 as nc

  # Identify closest TC
  dist_loc_cTC = [0.]*len(TCinfo["TClat"])
  for i in range(len(TCinfo["TClat"])):
    # Check if location is within distance to TC center
    dist_loc_cTC[i] = geodesic((lat,lon),(TCinfo["TClat"][i],
                                          TCinfo["TClon"][i])).km
  ind_closest_TC = dist_loc_cTC.index(min(dist_loc_cTC))
  dist_cTC = dist_loc_cTC[ind_closest_TC]
  TCradius = TCinfo["TCrad"][ind_closest_TC]

  # Check if within radius
  in_TC = 0
  TCname = ""
  if dist_cTC<TCradius:
    in_TC    = 1
    TCname   = str(nc.chartostring(
                  fTC.variables["name"][TCinfo["indTC"][
                  ind_closest_TC][0]]))

  return(dist_cTC,in_TC,TCname,TCradius)



#==================================================================
# Function to load land shape file
#==================================================================

def load_land(res='50m'):
  """
  Load land shape file

  Inputs:
   1) Optional resolution arguement, string format. Default = 50m.

  Outputs land shapoe file.

  Requires fiona 1.8.6 (conda install -c conda-forge fiona;
   https://pypi.org/project/Fiona/), cartopy 0.17.0
   (conda install -c conda-forge cartopy; 
   https://pypi.org/project/Cartopy/), and shapely 1.6.4
   (conda install -c conda-forge shapely; 
   https://pypi.org/project/Shapely/)
   
  """

  # Import libraries
  import fiona
  import cartopy.io.shapereader as shpreader
  import shapely.geometry as sgeom
  from shapely.prepared import prep
  
  # Pull in land shape files
  geoms = fiona.open(shpreader.natural_earth(resolution=res,
                           category='physical', name='land'))
  land_geom = sgeom.MultiPolygon([sgeom.shape(geom['geometry'])
                                for geom in geoms])
  land = prep(land_geom)
  
  # Check if point is over land and return information
  return(land)



#==================================================================
# Function to check if a point is over land
#==================================================================

def is_land(lon, lat, res='50m'):
  """
  Check if a location is over land. 

  Inputs:
   1,2) Longitude (x) and latitude (y) of position to check.
   3) Optional resolution arguement, string format. Default = 50m.

  Outputs True if over land. Outputs False if not.

  Requires fiona 1.8.6 (conda install -c conda-forge fiona;
   https://pypi.org/project/Fiona/), cartopy 0.17.0
   (conda install -c conda-forge cartopy; 
   https://pypi.org/project/Cartopy/), and shapely 1.6.4
   (conda install -c conda-forge shapely; 
   https://pypi.org/project/Shapely/)
   
  """

  # Import libraries
  import fiona
  import cartopy.io.shapereader as shpreader
  import shapely.geometry as sgeom
  from shapely.prepared import prep
  
  # Pull in land shape files
  geoms = fiona.open(shpreader.natural_earth(resolution=res,
                           category='physical', name='land'))
  land_geom = sgeom.MultiPolygon([sgeom.shape(geom['geometry'])
                                for geom in geoms])
  land = prep(land_geom)
  
  # Check if point is over land and return information
  return(land.contains(sgeom.Point(lon, lat)))



#==================================================================
# Calculate local solar time based on longitude
#==================================================================

def calc_local_sun_time(dt,lon):
  """
  Calculates the local solar time given a date and time and 
   location. Calculated as the UTC time plus an offset based 
   on the longitude. The offset is calculated by multiplying 
   the longitude by 24/360. 

  Inputs:
   1) The date and time (dt). dt is a string indicating UTC
    time in format YYYYMMDDhhmmss.
   2) The longitude position (lon). lon is the longitude 
    position as a float.

  Output is a string in YYYYMMDDhhmmss format indicating the
   local solar time.

  Note: This is not the actual local time. This should 
   typically only be used to calculate times for the 
   diurnal cycle.
  """

  # Import libraries
  import datetime

  # Make a datetime object for the current date and time
  timenow = datetime.datetime(int(dt[0:4]),int(dt[4:6]),
             int(dt[6:8]),int(dt[8:10]),int(dt[10:12]),0)

  # Calculate time with offset added
  newtime = timenow + datetime.timedelta(hours=lon*(24./360.))

  # Return time in same format as input
  return(str(newtime.year).zfill(4)+\
         str(newtime.month).zfill(2)+\
         str(newtime.day).zfill(2)+\
         str(newtime.hour).zfill(2)+\
         str(newtime.minute).zfill(2)+\
         str(newtime.second).zfill(2))



#==================================================================
# Create 2D dataframe given a list of coordinates and data
#==================================================================

def create_2d_dataframe(x,y,dx,dy,data):
  """
  Creates a dataframe given a set of pixel locations, pixel sizes,
   and the data in each pixel. Fills any coordinates that are in 
   the new dataframe but were not in the original lists with zeros.

  Inputs:
   1,2) Lists of pixel locations (indexes) x and y
   3,4) Scalars of pixel width (dx) and pixel height (dy)
   5) A corresponding list of the data values.

  Output is a 2D dataframe with x and y the indexes, and data as 
   the values

  Requires numpy 1.16.3 and pandas 0.24.2 (conda install -c 
   conda-forge pandas; https://pypi.org/project/pandas/)
  """

  # Import libraries
  import numpy as np
  import pandas as pd

  # Create coordinates
  ny = [round(i,2) for i in np.arange(min(y)-2*dy, 
                                      max(y)+2*dy,dy)]
  nx = [round(i,2) for i in np.arange(min(x)-2*dx, 
                                      max(x)+2*dx,dx)]

  # Generate dataframe
  df = pd.DataFrame(np.zeros([len(ny),len(nx)],float),
                                 [str(i) for i in ny],
                                 [str(i) for i in nx])

  # Populate dataframe
  for i in range(len(y)):
    df.loc[str(y[i]),str(x[i])] = data[i]

  # Return dataframe
  return(df)



#==================================================================
# Identify contiguous areas in 2d field 
#==================================================================

def label_wdiags(array):
  """
  This is a small adjustment to the scipy.ndimage.label() function
   to connect pixels at diagonals as the same feature. Original 
   scipy.ndimage.label() documentation can be found here: 
   https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.label.html

  Input:
   1) A 2D list, numpy array, or pandas dataframe with the 
    data you wish to identify features in.

  Output is a 2D array corresponding to data but with features 
   labelled as integers.

  Requires scipy 1.2.1 (conda install -c anaconda scipy; 
   https://pypi.org/project/scipy/)
  """

  # Import libraries
  import scipy.ndimage as ndimg

  # Setup a 2x2 binary array to account for diagonals
  s = ndimg.generate_binary_structure(2,2)

  # Run scipy.ndimage.label() with binary structure
  labels, numL = ndimg.label(array, structure=s)

  # Returns array with labels and the number of features
  return(labels,numL)



#==================================================================
# Find unique corner points from a list of pixel centers
#==================================================================

def find_corners(coords,dx,dy,dp=2):
  """
  Defines a list of coordinates for corners given a list of pixel
   centers.

  Inputs:
   1) A list of coordinate pairs i.e. [[x1,y1],[x2,y2],...]
   2,3) Scalars defining the width and height of the pixels
   4) A scalar defining the decimal places to round the resulting 
    coordinates. Optional, default = 2.

  Output is a list of coordinate pairs for the corners.

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/)

  Note, only returns unique corner coordinates (i.e. there are 
   no coordinates that are the same even if two pixels share a 
   corner).
  """

  # Import libraries
  import numpy as np

  # Define the distances to a corner from center
  dx2=dx/2; dy2=dy/2

  # Define all corners
  allcorners = [(round(x+dx2,dp),round(y+dy2,dp))
                              for x,y in coords] + \
               [(round(x-dx2,dp),round(y+dy2,dp)) 
                              for x,y in coords] + \
               [(round(x+dx2,dp),round(y-dy2,dp)) 
                              for x,y in coords] + \
               [(round(x-dx2,dp),round(y-dy2,dp)) 
                              for x,y in coords]

  # Get rid of all doubles and return list of corner coordinates
  return([list(rows) for rows in 
          np.unique(allcorners, axis=0)])



#==================================================================
# Find which points are inside a shape
#==================================================================

def points_in_shape(verts,points,widen=0):
  """
  Finds and returns all grid points within a list of vertices 
   defining a shape. 

  Input: 
   1) A list of tuples of the x,y coordinates of the vertices
   2) A list of tuples of the x,y coordinates of the points to 
       check whether they are in the shape
   3) The amount to widen the shape by. Most common use - due to
       computational issues, sometimes the shape doesn't include
       points on the edge. By widening the shape slightly, you 
       can ensure points are included. Positive values widen the 
       shape when the vertcies are ordered counterclockwise and
       negative widens the shape when the vertices are ordered
       clockwise. Default is 0 - no widening.

  Output 

  Requires maplotlib 3.0.3 (conda install -c conda-forge 
   matplotlib; https://pypi.org/project/matplotlib/)
  """

  # Import libraries
  from matplotlib.path import Path

  # Define the shape
  p = Path(verts)

  # Get list of booleans for if inside array
  grid = p.contains_points(points,radius=widen)

  # Return list of coordinate tuples for points inside shape
  return([points[i] for i, x in enumerate(grid) if x])

#==================================================================
# Find which points are inside a shape
#==================================================================

def calc_mjrmnrax(lons,lats):
  """
  Finds the major and minor axes of a set of points on earth

  Input: 
   1,2) Lists of longitude and latitude coordinates to fit the 
    major and minor axes too

  Output 

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/)
  """

  # Import libraries
  import numpy as np

  # Calculate center of mass of shape
  center = [np.mean(lons),np.mean(lats)]

  # Calculate eigen-value/vector pairs for largest piece
  eigvals, eigvecs = np.linalg.eig(np.cov((lons,lats)))

  # Calculate coordinates of axes
  lonseig = np.zeros((2,2)); latseig = np.zeros((2,2))
  lonseig[0,:],latseig[0,:] = np.vstack((
    center+(eigvals[0]*2)*eigvecs[0,:], 
    center-(eigvals[0]*2)*eigvecs[0,:])).T
  lonseig[1,:],latseig[1,:] = np.vstack((
    center+(eigvals[1]*2)*eigvecs[1,:], 
    center-(eigvals[1]*2)*eigvecs[1,:])).T

  # Calculate lengths and angles of the axes
  lengths = np.zeros(2); angles = np.zeros(2)
  lengths[0],angles[0] = calc_distandangle(
   lonseig[0,0],latseig[0,0],lonseig[0,1],latseig[0,1])
  lengths[1],angles[1] = calc_distandangle(
   lonseig[1,0],latseig[1,0],lonseig[1,1],latseig[1,1])

  # Calculate angle from north
  for i in range(len(angles)):
    if angles[i]>180:
      angles[i] = angles[i]-180

  # Calculate which is the major and minor axes
  mjrind1 = np.argmax(lengths)
  if hasattr(mjrind1, "__len__"): mjrind=mjrind1[0]
  else: mjrind=mjrind1
  if mjrind==0: mnrind=1
  if mjrind==1: mnrind=0

  # Assign axes
  mjrax_len = lengths[mjrind]
  mnrax_len = lengths[mnrind]
  mjrax_ang = angles[mjrind]
  mnrax_ang = angles[mnrind]

  # Return
  return(center,mjrax_len,mjrax_ang,mnrax_len,mnrax_ang)



#==================================================================
# Write netcdf variable
#==================================================================

def write_var(varname,long_name,description,dimname,dtype,units,
              fileout,datain,f):
  """
  Generic script for defining/opening a variable in a 
   netcdf file and then writing the associated data to the 
   variable."
  """

  try:
    datafile = fileout.createVariable(varname, dtype, (dimname))
  except:
    #print(varname + " already defined")
    datafile = fileout.variables[varname]
  #print("Writing "+varname+" to "+f)
  datafile.long_name   = long_name
  datafile.description = description
  datafile.units = units
  datafile[:] = datain 



#==================================================================
# Write netcdf groups as attribute and value pairs 
#==================================================================

def write_group(groupname,long_name,description,units,
                format1,fileout,datain,f):
   """
   Generic script for defining/opening a group in a netcdf
    file and then writing a python dictionary to that group
    as attribute:value pairs. Note: value can be a string, 
    integer, float, list, or anything else that can be 
    taken as an attribute in a netcdf file."
   """

   try:
     datafile  = fileout.createGroup(groupname)
   except:
     #print(groupname+" already defined")
     datafile = fileout.group[groupname]
   try:
     datadatafile = datafile.createGroup('data')
   except:
     #print("data"+groupname+" already defined")
     datadatafile = fileout.group[groupname].groups["data"]

   #print("Writing "+groupname+" to "+f)
   datafile.long_name   = long_name
   datafile.description = description
   datafile.units = units
   datafile.format = format1
   for k,v in datain.items():
     setattr(datadatafile, k,  v)



#==================================================================
# End
#==================================================================
