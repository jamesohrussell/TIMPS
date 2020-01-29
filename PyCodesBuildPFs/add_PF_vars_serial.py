#==========================================================
# Add variables to IMERG PF netcdf files.
#
# Must provide in namelist:
#  * Directory and filename identifiers for PF files
#  * Directory and filename identifiers for IBtracs data
#     (only if TC information is desired)
#  * Variables desired
#  * Grid spacing of input data pixels in degrees 
#     (IMERG = 0.1, only required if area or volumetric 
#      rain rate is desired)
#  * Subset information if desired
#
# James Russell 2019
#==========================================================

# Import python libraries (do not change)
import numpy as np
from netCDF4 import Dataset
import glob
import datetime
import driver_addvars as da
from joblib import Parallel, delayed
import PF_functions as PFfunc
import time as tm
import os
import fiona
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep

#==========================================================
# Namelist
#==========================================================

# Directory and filename for PF files
datadir = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_100719/IMERGPFs/"
fileid  = "IPF_"

# Subset (for certain range of dates) 
ssdat = False
date1 = "20180830"
date2 = "20180902"
ssobs = True
obid1 = "118772"
obid2 = "118773"

# Variables desired
addmaxrainrate    = True # Maximum rain rate
addmeanrainrate   = True # Mean rain rate
addmedianrainrate = True # Median rain rate
addstddevrainrate = True # Standard deviation of the rain rates
addarea           = True # Area of the PF
addvolrainrate    = True # Volumetric rain rate
addpropagation    = True # Propagation characteristics
addtimevars       = True # Add time variables
addTCinfo         = True # Flags indicating proximity to TC center
addlandinfo       = True # Flags indicating locations over land

# Inputs for specific variables
# Grid spacing in degrees lon, lat (only required for area/volrainrate)
dx = 0.1
dy = 0.1
# Directory and filename of TC data (only required for TC information)
dataTCdir = "/uufs/chpc.utah.edu/common/home/varble-group2/IBTrACS/"
fileTCid  = "IBTrACS.ALL.v04r00.nc"

#==========================================================
# Initialize timer
#==========================================================

start = tm.time()

#==========================================================
# Generate list of IPF files to process
#==========================================================

# Reads directory and file names
print("Generating file list")
if ssdat:

  # Generate a list of filenames with dates to search for
  start = datetime.datetime.strptime(date1,"%Y%m%d")
  end = datetime.datetime.strptime(date2,"%Y%m%d")
  datearr = (start + datetime.timedelta(days=x) for x in range(0,(end-start).days))
  filen = []
  for dateobj in datearr:
    filen.append(datadir+fileid+"*"+dateobj.strftime("%Y%m%d")+"*")

  # Reads directory and filenames
  filenamesrun = []
  for f in filen:
      filenamesrun = filenamesrun+glob.glob(f)

elif ssobs:

  obids = [i for i in range(int(obid1),int(obid2))]
  filen = []
  for i in range(len(obids)):
    filen.append(datadir+fileid+str(obids[i])+"*")

  # Reads directory and filenames
  filenamesrun = []
  for f in filen:
      filenamesrun = filenamesrun+glob.glob(f)

else:
  # Find all files in directory
  filenamesrun = glob.glob(datadir+fileid+'*')

#==========================================================
# Loop over IPF files in parrallel
#==========================================================

# Loop over PFs
print("Begin loop over IPFs")

for f in filenamesrun:
  print("Working with file "+str(f))

  # Open netcdf file and assign data
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

  # Preallocate arrays
  if addmaxrainrate:
    maxrainrate1  = [0.]*len(datakeys)

  if addmeanrainrate:
    meanrainrate1 = [0.]*len(datakeys)

  if addmedianrainrate:
    medianrainrate1 = [0.]*len(datakeys)

  if addstddevrainrate:
    stddevrainrate1 = [0.]*len(datakeys)

  if addarea:
    area1 = [0.]*len(datakeys)

  if addvolrainrate:
    volrainrate1 = [0.]*len(datakeys)

  if addpropagation:
    propspd1 = [0.]*len(datakeys)
    propdir1 = [0.]*len(datakeys)

  if addtimevars:
    # Calculate normalized time variable
    np.seterr(divide='ignore', invalid='ignore')
    normalizedtime1 = [i for i in range(len(datakeys))]/np.float32(len(datakeys)-1)

  if addTCinfo:
    dist_cPF_cTC1 = [0.]*len(datakeys)
    cPF_in_TC1    = [0]*len(datakeys)
    TCname_cPF1   = [""]*len(datakeys) 
    TCrad_cPF1    = [0.]*len(datakeys)
    lPF_in_TC1    = {}
    dist_lPF_cTC1 = {}
    TCname_lPF1   = {}
    TCrad_lPF1    = {}
    writeTCdata   = False

  if addlandinfo:
    cPF_over_land = [0]*len(datakeys)
    lPF_over_land = {}

#==========================================================
# Begin loop over times
#==========================================================

  c = 0
  for k in datakeys:
    print(str(c)+": Getting data and calculating new variables for "+k)

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

#==========================================================
# Add land surface information
#=========================================================

    if addlandinfo:

      print("Calculating land information")

      # Pull in land shape files
      geoms = fiona.open(shpreader.natural_earth(resolution='50m',
                           category='physical', name='land'))
      land_geom = sgeom.MultiPolygon([sgeom.shape(geom['geometry'])
                                for geom in geoms])
      land = prep(land_geom)

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

#==========================================================
# Calculate simple variables to add
#==========================================================

    if addmaxrainrate:
      print("Calculating max rain rate")
      maxrainrate1[c] = np.amax(instrainnzk)

    if addmeanrainrate:
      # Calculate mean rain rates but exclude pixels with no rain
      print("Calculating mean rain rate")
      meanrainrate1[c]   = np.mean(instrainnzk)

    if addmedianrainrate:
      # Calculate median rain rates but exclude pixels with no rain
      print("Calculating median rain rate")
      medianrainrate1[c] = np.median(instrainnzk)

    if addstddevrainrate:
      # Calculate standard deviation of rain rates but exclude pixels with no rain
      print("Calculating standard deviation of rain rates")
      stddevrainrate1[c] = np.std(instrainnzk)

#==========================================================
# Add area or volumetric rain rate information
#==========================================================

    if (addarea) or (addvolrainrate):
      # Calculate area
      print("Calculating area and volumetric rain rate")

      area1[c],volrainrate1[c] = PFfunc.calc_area_and_volrainrate(
                      lonsnzk,latsnzk,instrainnzk,float(dx),float(dy))

#==========================================================
# Add propagation
#==========================================================

    if addpropagation:
      # Calculate speed and direction of motion of PF centroid
      print("Calculating propagation information")

      # Account for objects that aren't more than one time
      if len(lats)<2:
        propspd1[c] = 0
        propdir1[c] = 0

      else:

        # If only one time in PF, set as 0
        #print("Calculating propagation speed and direction")
        # Forward difference for first time
        if c==0:

          propspd1[c],propdir1[c] = PFfunc.calc_propagation(
                                   datadtim[c],dataclat[c],dataclon[c],
                                   datadtim[c+1],dataclat[c+1],dataclon[c+1])

        # Backward difference for last time
        elif c==len(dataclon)-1:
          propspd1[c],propdir1[c] = PFfunc.calc_propagation(
                                   datadtim[c-1],dataclat[c-1],dataclon[c-1],
                                   datadtim[c],dataclat[c],dataclon[c])

        # Centered difference for all other times
        else:
          propspd1[c],propdir1[c] = PFfunc.calc_propagation(
                                   datadtim[c-1],dataclat[c-1],dataclon[c-1],
                                   datadtim[c+1],dataclat[c+1],dataclon[c+1])
     
#==========================================================
# Add TC information
#==========================================================

    if addTCinfo:
      print("Calculating tropical cyclone information")
      
      # Read tropical cyclone data
      fTC = Dataset(dataTCdir+fileTCid) 

      # Interpolate times and latitudes and get radius
      TCinfo = PFfunc.interp_TC(datadtim[c],fTC)

      if TCinfo["TC"]:

        # Check if center is close to TC
        dist_cPF_cTC1[c],cPF_in_TC1[c],TCname_cPF1[c],TCrad_cPF1[c] \
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
        lPF_in_TC1[k]    = TCbinnow
        dist_lPF_cTC1[k] = TCdistnow
        TCname_lPF1[k]   = TCnamenow
        TCrad_lPF1[k]    = TCradnow

        # Set boolean to indicate if writing data is required
        if (cPF_in_TC1[c]==1) or (1 in lPF_in_TC1[k]):
          writeTCdata = True

      else:
        #print("No tropical cyclones at this time")
        writeTCdata = False

#==========================================================
# End loops over objects and times
#==========================================================
      
    # Advance counter
    c = c + 1

  fd.close()

#==========================================================
# Write all data to file
#==========================================================

  print("Writing data")

  fileout = Dataset(f,'a')

  # Write time variables to file
  if addtimevars:
    description = "A fractional time with 0 the first time of the PF and 1 the last time of the PF"
    PFfunc.write_var("normalizedtime","Normalized time",description,
                 "time",np.float64,"",fileout,normalizedtime1,f)

  # Write max rain rate to file
  if addmaxrainrate:
    description = "Maximum rain rate within PF"
    PFfunc.write_var("maxrainrate","Max rain rate",description,
                 "time",np.float64,"mm/hr",fileout,maxrainrate1,f)

  # Write mean rain rate to file
  if addmeanrainrate:
    description = "Mean rain rate within PF excluding pixels with zero rain rate"
    PFfunc.write_var("meanrainrate","Mean rain rate",description,
                 "time",np.float64,"mm/hr",fileout,meanrainrate1,f)

  # Write median rain rate to file
  if addmedianrainrate:
    description = "Median rain rate within PF excluding pixels with zero rain rate"
    PFfunc.write_var("medianrainrate","Median rain rate",description,
                 "time",np.float64,"mm/hr",fileout,medianrainrate1,f)

  # Write standard deviation rain rate to file
  if addstddevrainrate:
    description = "Standard deviation of rain rate within PF excluding pixels with zero rain rate"
    PFfunc.write_var("stddevrainrate","Standard deviation of rain rate",
                 description,"time",np.float64,"mm/hr",fileout,
                 stddevrainrate1,f)

  # Write area to file
  if addarea:
    description = "Area within PF excluding pixels with zero rain rate"
    PFfunc.write_var("area","Area",description,"time",np.float64,
                 "km^2",fileout,area1,f)

  # Write volumetric rain rate to file
  if addvolrainrate:
    description = "Volumetric rain rate within PF excluding pixels with zero rain rate. Calculated as the sum over all pixels (with non-zero rain rate) of the area of the pixel multipled by the rain rate of the pixel."
    PFfunc.write_var("volrainrate","Volumetric rain rate",description,"time",
                 np.float64,"mm hr^-1 km^2",fileout,volrainrate1,f)

  # Write propagation to file
  if addpropagation:
    description = "Calculated as the geodesic distance travelled by centroid divided by time taken"
    PFfunc.write_var("propspd","Propagation speed",description,"time",
                 np.float64,"m/s",fileout,propspd1,f)

    description = "Calculated as the clockwise angle from north the centroid is moving toward."
    PFfunc.write_var("propdir","Propagation direction",description,"time",
                 np.float64,"degrees",fileout,propdir1,f)

  # Write land information
  if addlandinfo:
    description = "1 if center of PF is over land. 0 if not."
    PFfunc.write_var("cPF_over_land","Center of PF over land",description,"time",
                 int,"",fileout,cPF_over_land,f)

    format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."
    description = "Binary indicating if a PF location is over land (1 = over land, 0 = not over land)"
    PFfunc.write_group("lPF_over_land","Location over land",description,"",
                     format1,fileout,lPF_over_land,f)

  # Write TC information to file
  if addTCinfo:
    
    if writeTCdata: 

    # If any part of PF IS within TC radius set TC attribute as True and write
    #  other variables describing which parts of PF are within TC radius 
      #print("PF within TC, writing data")
      fileout.within_TC = "True"

      # Write cPF_in_TC
      description = "If PF center is within a maximum radius of TC center this is 1, if not, this is 0."
      PFfunc.write_var("cPF_in_TC","Proximity of PF center to TC",description,
                 "time",int,"",fileout,cPF_in_TC1,f)

      # Write dist_cPF_cTC
      description = "Calculated as the geodesic distance of the PF center to TC center"
      PFfunc.write_var("dist_cPF_cTC","Distance of PF center to TC center", 
                   description,"time",np.float64,"km",fileout,dist_cPF_cTC1,f)

      # Write TCradius
      description = "Calculated by finding the largest radius of a TC from IBtracs and doubling."
      PFfunc.write_var("TCrad_cPF","Radius of TC",description,"time",
                   np.float64,"km",fileout,TCrad_cPF1,f)

      # Write list of TCname_cPF as an attribute of PF file
      fileout.TCname_cPF = TCname_cPF1
     
      # Create or open a group for lPF_in_TC
      format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."
      description = "Binary indicating if a PF location is within a TC (1 = within TCradius, 0 = not within TC radius)"
      PFfunc.write_group("lPF_in_TC","Location within TC",description,"",
                     format1,fileout,lPF_in_TC1,f)

      # Create or open a group for dist_lPF_cTC
      description = "Calculated as geodesic distance to TC center from PF location"
      PFfunc.write_group("dist_lPF_cTC","Distance to TC center from PF location",
                     description,"km",format1,fileout,dist_lPF_cTC1,f)

      # Create or open a group for TCname_lPF
      description = "Name of TC if given location is within TC"
      Pfunc.write_group("TCname_lPF","Name of TC if given location is within TC",
                     description,"",format1,fileout,TCname_lPF1,f)

      # Create or open a group for TCrad_lPF
      description = "Calculated by finding the largest radius of a TC from IBtracs and doubling."
      PFfunc.write_group("TCrad_lPF","Radius of TC",
                     description,"",format1,fileout,TCrad_lPF1,f)

    else:
    # If not, just set TC attribute as False
      #print("PF not within TC")
      fileout.within_TC = "False"

  # Close current file
  fileout.close()

#==========================================================
# Final clean up tasks
#==========================================================

# End timer
end = tm.time()
print("Program took "+str(end - start)+"s")

#==========================================================
# End code
#==========================================================
