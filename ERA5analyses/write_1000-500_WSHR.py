
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
lon1    = -180
lon2    = 179.75

#==================================================================
# Import libraries
#==================================================================

import ERA5_functions as E5fns
from netCDF4 import Dataset
import numpy as np

#==================================================================
# Get variable
#==================================================================

SPHU,times,levs,lats,lons = E5fns.get_E5_ss_4D(datadir,fileid,
 varname,timestr1,timestr2,lev1,lev2,lat1,lat2,lon1,lon2)

#==================================================================
# Average across levels
#==================================================================

Q = np.mean(SPHU,axis=1)

print(np.shape(Q))

#==================================================================
# Write to file
#==================================================================












