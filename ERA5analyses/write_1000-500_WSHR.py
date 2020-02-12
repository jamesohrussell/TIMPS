import ERA5_functions as E5fns
from netCDF4 import Dataset

datadir = "/uufs/chpc.utah.edu/common/home/varble-group1/ERA5/"
fileid  = "ERA5.VWND.20"
varname = "V"
timestr1= "2018-06-30 11:00:00"
timestr2= "2018-07-02 15:45:00"
lat1    = -10
lat2    = 10
lon1    = -20
lon2    = 30

files,timi1,timi2,times = E5fns.get_E5_ss_4D_fiti(datadir,fileid,timestr1,timestr2)

fh = Dataset(files[0])

loni,lati,lons,lats = E5fns.get_E5_ss_4D_coords(fh,lon1,lon2,lat1,lat2)

print(loni)
print(lati)
print(lons)
print(lats)
