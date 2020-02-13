
#==================================================================
# Namelist
#==================================================================

datadir = "/uufs/chpc.utah.edu/common/home/varble-group1/ERA5/"
fileid  = "ERA5.SPHU."
varname = "Q"
# Times must be in format "YYYY-MM-DD hh"
timestr1= "2018-06-01 00"
timestr2= "2018-06-30 23"
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

import ERA5_functions as E5fns
from netCDF4 import Dataset
import numpy as np
import re
import time as tm
import TIPS_functions as fns

#==================================================================
# Get variable
#==================================================================

SPHU,times,levs,lats,lons = E5fns.get_E5_ss_4D(datadir,fileid,
 varname,timestr1,timestr2,lev1,lev2,lat1,lat2,lon1,lon2)

#==================================================================
# Average across levels
#==================================================================

Q1 = np.mean(SPHU,axis=1)

#============================================================
# Write NetCDF file for object
#============================================================

date1 = re.findall(r"[\w']+", timestr1)
date2 = re.findall(r"[\w']+", timestr2)

# Creating netcdf file to write to
filename = "ERA5.SPHU_1000-850hPamean."+\
      date1[0]+date1[1]+date1[2]+date1[3]+"_"+\
      date2[0]+date2[1]+date2[2]+date2[3]+".nc"

print("Creating file: "+datadir+filename)
fileout = Dataset(datadir+filename,"w",format="NETCDF4")

# Global Attributes
fileout.long_name   = "Specific humidity (1000-850 hPa mean)"
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
fns.write_var("Q","Specific humidity",
 "Specific humidity averaged between 1000-850 hPa",
 ("time","latitude","longitude"),np.float64,"",fileout,
 Q1,datadir+filename)

#============================================================
# End
#============================================================
