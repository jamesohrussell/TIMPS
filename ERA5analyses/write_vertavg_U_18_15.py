
#==================================================================
# Namelist
#==================================================================

datadir = "/uufs/chpc.utah.edu/common/home/varble-group1/ERA5/winds/"
fileid  = "ERA5.UWND."
varname = "U"
varfile = "UWND"
varfileln = "Zonal wind"
varfileds = "Zonal wind (1000-850 hPa mean)"
# Times must be in format "YYYY-MM-DD hh"
timestr1= ["2015-01-01 00","2015-02-01 00","2015-03-01 00",
           "2015-04-01 00","2015-05-01 00","2015-06-01 00",
           "2015-07-01 00","2015-08-01 00","2015-09-01 00",
           "2015-10-01 00","2015-11-01 00","2015-12-01 00"]
timestr2= ["2015-01-31 23","2015-02-28 23","2015-03-31 23",
           "2015-04-30 23","2015-05-31 23","2015-06-30 23",
           "2015-07-31 23","2015-08-31 23","2015-09-30 23",
           "2015-10-31 23","2015-11-30 23","2015-12-31 23"]
# Levels are 1-1000 (hPa)
lev1    = 850
lev2    = 1000
# Latitudes are in -90-90
lat1    = -90
lat2    = 90
# Longitudes can be in 0-359.75 or -180-179.75
lon1    = 0
lon2    = 359.75

tunits   = "hours since 1900-01-01 00:00:00"
qunits   = "kg kg**-1"

#==================================================================
# Import libraries
#==================================================================

import sys
sys.path.insert(0,'/uufs/chpc.utah.edu/common/home/u0816744/TIPS_git/PyCodesBuildTIPS/')
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
  levi,levs = E5fns.get_E5_ss_4D_levels(fh,lev1,lev2)

  # Get the variable 
  var = E5fns.get_E5_ss_4D_levavgvar(
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
  fns.write_var(varname,varfileln,varfileds,
   ("time","latitude","longitude"),np.float64,qunits,fileout,
   var,datadir+filename)

#==================================================================
# End
#==================================================================
