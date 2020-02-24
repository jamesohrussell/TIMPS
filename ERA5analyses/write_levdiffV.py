
#==================================================================
# Namelist
#==================================================================

datadir = "/uufs/chpc.utah.edu/common/home/varble-group1/ERA5/"
fileid  = "ERA5.VWND."
varname = "V"
varfile = "VSHR"
varfileln = "Meridional Wind Shear"
varfileds = "Meridional Wind Shear (200 - 850 hPa)"
# Times must be in format "YYYY-MM-DD hh"
timestr1= ["2017-01-01 00","2017-02-01 00","2017-03-01 00",
           "2017-04-01 00","2017-05-01 00","2017-06-01 00",
           "2017-07-01 00","2017-08-01 00","2017-09-01 00",
           "2017-10-01 00","2017-11-01 00","2017-12-01 00"]
timestr2= ["2017-01-31 23","2017-02-28 23","2017-03-31 23",
           "2017-04-30 23","2017-05-31 23","2017-06-30 23",
           "2017-07-31 23","2017-08-31 23","2017-09-30 23",
           "2017-10-31 23","2017-11-30 23","2017-12-31 23"]
# Levels are 1-1000 (hPa)
lev1    = 200
lev2    = 850
# Latitudes are in -90-90
lat1    = -90
lat2    = 90
# Longitudes can be in 0-359.75 or -180-179.75
lon1    = 0
lon2    = 359.75

tunits   = "hours since 1900-01-01 00:00:00"
dunits   = "ms^-1"

#==================================================================
# Import libraries
#==================================================================

import ERA5_functions as E5fns
from netCDF4 import Dataset
import numpy as np
import re
import time as tm
import TIPS_functions as fns

#==================================================================
# Get variable
#==================================================================

for i in range(0,len(timestr1)):

  print("Working on dates "+timestr1[i]+" -- "+timestr2[i])

  # Get files, times, and time indices for first and last file
  files,timi,times = E5fns.get_E5_ss_4D_fiti(datadir,fileid,
                                        timestr1[i],timestr2[i])

  # Read first file and get lats, lons, and their indices 
  fh = Dataset(files[0])
  lati,loni,lats,lons = E5fns.get_E5_ss_4D_coords(fh,lat1,lat2,
   lon1,lon2)

  # Get levels and their indices 
  levi = E5fns.get_E5_ss_4D_2levels(fh,lev1,lev2)

  # Get the variable 
  var = E5fns.get_E5_ss_4D_levdiffvar(
   files,varname,timi,levi,lati,loni)

#==================================================================
# Write NetCDF file for object
#==================================================================

  date1 = re.findall(r"[\w']+", timestr1[i])
  date2 = re.findall(r"[\w']+", timestr2[i])

  # Creating netcdf file to write to
  filename = "ERA5."+varfile+"_"+\
        str(lev2)+"-"+str(lev1)+"hPamean."+\
        date1[0]+date1[1]+date1[2]+date1[3]+"_"+\
        date2[0]+date2[1]+date2[2]+date2[3]+".nc"

  print("Creating file: "+datadir+filename)
  fileout = Dataset(datadir+filename,"w",format="NETCDF4")

  # Global Attributes
  fileout.long_name   = varfileln
  fileout.description = ""
  fileout.history     = 'Created ' + tm.ctime(tm.time())
  fileout.source      = 'ERA5: NCAR RDA'

  # Create dimensions in file
  time = fileout.createDimension('time',len(times))
  lat  = fileout.createDimension('latitude',len(lats))
  lon  = fileout.createDimension('longitude',len(lons))

  # Create dimension variables in file
  fns.write_var("time","Time","","time",np.float64,
                 tunits,fileout,times,datadir+filename)

  fns.write_var("latitude","Latitude","","latitude",np.float64,
                "degreesN",fileout,lats,datadir+filename)

  fns.write_var("longitude","Longitude","","longitude",np.float64,
                "degreesE",fileout,lons,datadir+filename)

  # Create variables in file
  fns.write_var(varfile,varfileln,varfileds,
   ("time","latitude","longitude"),np.float64,dunits,fileout,
   var,datadir+filename)

#==================================================================
# End
#==================================================================
