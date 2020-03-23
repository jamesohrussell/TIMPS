#==========================================================
# Script reads IPF data in a folder and plots all tracks on
#   a map of the tracked domain.
#
# James Russell 2019
#==========================================================

# Import python libraries (do not change)
import numpy as np
from netCDF4 import Dataset
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.pyplot as plt
import cmaps
from mpl_toolkits.axes_grid1 import make_axes_locatable

# For plotting
plot  = True # Turn on or off plot (mostly for debugging)
labnx = 8   # Number of x-labels between/incl. lonW,lonE
labny = 9    # Number of y-labels between/incl. latS,latN
dx    = 0.5  # Grid-spacing in x (0.1 for IMERG)
dy    = 0.5  # Grid spacing in y (0.1 for IMERG)
latmn = -9   # Min lat for plotting
latmx = 29   # Max lat for plotting
lonmn = -69  # Min lon for plotting
lonmx = -1  # Max lon for plotting
mincl = 0    # Min for colorbar
maxcl = 0.0075   # Max for colorbar

#==========================================================
# Read data
#==========================================================

# Open file
dataset0 = Dataset("ugdatabytime_map0.5_2014-2018.nc")
# Read variables
avga = dataset0.variables["bins"][:,:]

lat = np.arange(-10,30.+dy,dy)
lon = np.arange(-70,0+dx,dx)

latN = np.max(lat)
latS = np.min(lat)
lonE = np.max(lon)
lonW = np.min(lon)

#==========================================================
# Make plot
#==========================================================

print("Plotting data")

# Generate center points for pcolor
lons,lats = np.meshgrid(lon,lat)

# Make figure
ax = plt.figure(num=None, figsize=(8, 6), dpi=200, facecolor='w', edgecolor='k')
#ax = plt.subplot()
  
# set up a map
extent = [lonW, lonE, latS, latN]
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent(extent,crs=ccrs.PlateCarree())
ax.coastlines(resolution='50m')

# Make plot
cs = plt.pcolor(lons,lats,avga,\
                    transform=ccrs.PlateCarree(),\
                    cmap=cmaps.WhiteYellowOrangeRed)#,\
#                    vmin=mincl,vmax=maxcl)

# Format tick marks
ax.set_xticks(np.linspace(lonW,lonE,num=labnx), crs=ccrs.PlateCarree())
ax.set_yticks(np.linspace(latS,latN,num=labny), crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

cbar = plt.colorbar(cs, ax=ax, orientation="vertical", pad=0.01, fraction=0.027)   
cbar.set_label("# of occurences of upscale growth",rotation=270,labelpad=10)

plt.xlim(lonmn, lonmx)
plt.ylim(latmn, latmx)
#plt.title("Upscale growth for Jun-Sep 2018 ")
plt.savefig('ugbytime_14-18.png')

