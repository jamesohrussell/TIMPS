#==========================================================
# Script reads IMERG data files and tracked FiT objects,
# then reworks data such that all are on the same grid,
# and plots them with IMERG data pixels and black outlines
# denoting FiT object(s) (or thresholds for testing).
#
# In development:
#  * Plot tracks of IPF centroids.
#
# Must provide in namelist:
#  * Directory and filename identifiers for IMERG files 
#      and FiT output netcdf files (the latter has no 
#      effect if only plotting IMERG)
#  * Domain bounds for data region if not global (needs to
#      be exactly the same as region tracked using FiT)
#  * Gaussian smoothing information for IMERG (should be 
#      the same as for FiT input files)
#  * Various plotting options
#
# James Russell 2019
#==========================================================

# Import python libraries (do not change)
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import glob
import time
import scipy.ndimage as ndimage
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps
import h5py
from collections import Counter
import csv
import datetime

#==========================================================
# Namelist
#==========================================================

# Directory and filename for IMERG data (specify format --- only hdf5 or netcdf4 files)
# Should obtain from create_FiT_input_files.py
datadirIM = "/uufs/chpc.utah.edu/common/home/varble-group2/IMERG/"
fileidIM  = "3B-HHR.MS.MRG.3IMERG."
datahdf5  = True
datanc4   = False

# Date and time range
starttime = "20180601"
endtime   = "20181001" # Actual last day is day before

# Directory and filename for FiT output data files
datadirFi = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_100719/out_files/"
fileidFi1 = "IMERG_tracked_"
fileidFi2 = "_4Dobjects.nc"
fileidFi3 = "IMERG_tracked_4Dobject_tree.txt"

# Directory and filename for IPFs
#datadirPF = "/uufs/chpc.utah.edu/common/home/varble-group2/james/IMERGPrecipitationFeatures/"
#fileidPF  = "IPF_"

# Subset regions (if global set ssreg=False)
# Should obtain from create_FiT_input_files.py
ssreg  = True
latN   = 20
latS   = -5
lonW   = -60
lonE   = -15
ssname = "CapVer"

# Gaussian smoothing (No smoothing if False)
# Should obtain from create_FiT_input_files.py
smooth = False
gaussian = False
stddev   = 1.5  # Standard deviation for gaussian smoothing
uniform  = False
width    = 5    # Number of points to spread running average

# if plotobj=True, plots object outlines 
# if idobj=False, outlines all objects after obid1
# if idobj=True, outlines only object with obid1
# if longobj=True, plots nlong longest lived objects
plotobj = True
idobj   = False
obid1   = 100039
longobj = False
nlong   = 10

# If plotthr=True plots thresholds used for FiT
# This option should be False if plotobj = True as this 
#   will just overwrite the contours there. 
# Should obtain thresholds from create_FiT_input_files.py
plotthr = False
tholds  = [0.8,3.0,9.0,27.0]

# For plotting
plot  = True # Turn on or off plot (mostly for debugging)
labnx = 10    # Number of x-labels between/incl. lonW,lonE
labny = 6    # Number of y-labels between/incl. latS,latN
dx    = 0.1  # Grid-spacing in x (0.1 for IMERG)
dy    = 0.1  # Grid spacing in y (0.1 for IMERG)
latmn = -5   # Min lat for plotting
latmx = 20   # Max lat for plotting
lonmn = -60 # Min lon for plotting
lonmx = -15  # Max lon for plotting
mincl = -0.5 # Min for colorbar
maxcl = 16.5 # Max for colorbar
#tracks= True # Plot tracks?

#==========================================================
# Function to read text file
#==========================================================

def read_FiT_txt(dirandfile):
  "This function opens the FiT .txt output file and reads in the data"

  # Open file and give it the handle f
  f = open(dirandfile,'r')

  # Use numpy function to read data from .eol file
  data = np.genfromtxt(f, dtype="str",delimiter="\t",skip_header=1)

  return(data)

  # Close file
  f.close()

#==========================================================
# Begin loop
#==========================================================

# Generate a list of filenames with dates to search for
print("Finding all files in date range")
start = datetime.datetime.strptime(starttime, "%Y%m%d")
end = datetime.datetime.strptime(endtime, "%Y%m%d")
datearr = (start + datetime.timedelta(days=x) for x in range(0, (end-start).days))
filen = []
for dateobj in datearr:
  filen.append(datadirIM+fileidIM+dateobj.strftime("%Y%m%d")+"*")

# Reads directory and filenames
n = 0
for f in filen:
  if n==0:
    filenames = glob.glob(f)
  else:
    filenames = filenames+glob.glob(f)
  n=n+1
lenIdir    = len(datadirIM)

# Begin loop over all times
for i in range(len(filenames)):

  print("Reading data from "+filenames[i])

#==========================================================
# Read IMERG data (NetCDF4)
#==========================================================

  if datanc4:

    # Open file
    datasetIM = Dataset(filenames[i])

    # Read variables
    rain = datasetIM.variables["precipitationCal"][:,:,:]

#==========================================================
# Read IMERG data (hdf5)
#==========================================================

  if datahdf5:

    # Open file
    datasetIM = h5py.File(filenames[i],'r')

    # Read variables
    rain = datasetIM["/Grid/precipitationCal"][:,:,:]

#==========================================================
# Process IMERG data
#==========================================================

  # Change dimension order of variable
  raint = np.transpose(rain, (0, 2, 1))
  del(rain)
  rain = raint
  del(raint)

#==========================================================
# Read FiT tracking data
#==========================================================

  if plotobj:
    # Open file
    datasetFi = Dataset(datadirFi+fileidFi1+str(i).zfill(5)+fileidFi2)
  
    # Read variables
    obid = np.squeeze(datasetFi.variables["value"][:,:,:])

#==========================================================
# Process FiT tracking data
#==========================================================

    if idobj:
      # Make obbin only 1 for specific object id
      obbin = np.where(obid==obid1,1,0)
    elif longobj:
      # Select largest nlngst objects
      datatxt = read_FiT_txt(datadirFi+fileidFi3)
      obidsdrs = Counter([int(i) for i in datatxt[:,0]])
      obidssrt = sorted(obidsdrs, key=obidsdrs.get,reverse=True)
      obids = obidssrt[0:nlong]

#==========================================================
# Take time from IMERG file name
#==========================================================
  
  ye = str(filenames[i][lenIdir+21:lenIdir+25])
  mo = str(filenames[i][lenIdir+25:lenIdir+27])
  da = str(filenames[i][lenIdir+27:lenIdir+29])
  ho = str(filenames[i][lenIdir+31:lenIdir+33])
  mi = str(filenames[i][lenIdir+33:lenIdir+35])

#==========================================================
# Subset data
#==========================================================

  if ssreg:

    print("Subsetting variable by area")
    
    # Read in coordinate variables
    if datanc4:
      lat = datasetIM.variables["lat"][:]
      lon = datasetIM.variables["lon"][:]

    if datahdf5:
      lat = datasetIM["Grid/lat"][:]
      lon = datasetIM["Grid/lon"][:]

    # Find indices of closest coordinates
    idyS = (np.abs(lat - latS)).argmin()
    idyN = (np.abs(lat - latN)).argmin()
    idxW = (np.abs(lon - lonW)).argmin()
    idxE = (np.abs(lon - lonE)).argmin()

    # Subet main dataset by indices
    latsub  = lat[idyS:idyN]
    lonsub  = lon[idxW:idxE]
    rainsub = rain[:,idyS:idyN,idxW:idxE]

    # Replace old vars with new subsetting vars
    del(rain)
    del(lat)
    del(lon)

    rain = rainsub
    lat = latsub
    lon = lonsub

    del(rainsub)
    del(latsub)
    del(lonsub)

#==========================================================
# Smooth data
#==========================================================

  if smooth:
    print("Smoothing rainfall")
    if gaussian:
      rain = ndimage.gaussian_filter(rain, sigma=stddev)
    elif uniform:
      rain = ndimage.uniform_filter(rain,size=width)

#==========================================================
# Make plot
#==========================================================

  if plot:
    print("Plotting data")

    # Generate center points for contourf
    lonc,latc = np.meshgrid(lon,lat) 

    # Generate point bounds for pcolorf
    lats, lons = np.mgrid[slice(latS, latN + dy, dy),
                    slice(lonW, lonE + dx, dx)]

    # Make figure
    ax = plt.figure(num=None, figsize=(8, 6), dpi=200, facecolor='w', edgecolor='k')
    #ax = plt.subplot()
  
    # set up a map
    extent = [lonW, lonE, latS, latN]
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(extent,crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m')

    # Format tick marks
    ax.set_xticks(np.linspace(lonW,lonE,num=labnx), crs=ccrs.PlateCarree())
    ax.set_yticks(np.linspace(latS,latN,num=labny), crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    # Make plot
    cs = plt.pcolormesh(lons,lats,np.squeeze(rain),\
                        transform=ccrs.PlateCarree(),\
                        cmap=cmaps.prcp_1,\
                        vmin=mincl,vmax=maxcl)
    cbar = plt.colorbar(cs, ax=ax, orientation="horizontal", pad=.1, fraction=0.06)   
    cbar.set_label("IMERG Precipitation (mm/hr)")

    # Contour outline of objects
    if plotobj:

      cntrclr = ['r','g','b','m','y','#1f77b4','#ff7f0e','#9467bd','#8c564b','#7f7f7f']        
     
      if longobj:
        for p in range(0,nlong):
          obbin = np.where(obid==obids[p],1,0)
          plt.contour(lonc, latc, np.squeeze(obbin), [1], colors=cntrclr[p])
          del(obbin)
      else:
        obidsnow = np.unique(obid)
        obidsnow = [int(i) for i in obidsnow]
        obidsnow = list(filter((0.).__ne__, obidsnow))
        for i in obidsnow:
          obbinnow = np.where(obid==i,1,0)
          plt.contour(lonc, latc, np.squeeze(obbinnow), [1], colors=cntrclr[i%10])

    # Contour outline of thresholds
    if plotthr:
      plt.contour(lonc, latc, np.squeeze(rain), tholds, colors='k')

    plt.xlim(lonmn, lonmx)
    plt.ylim(latmn, latmx)
    plt.annotate(ye+"/"+mo+"/"+da+" "+ho+":"+mi, xy=(0.78, 0.92), xycoords='axes fraction')
    plt.savefig('figs/IMERG_'+ssname+"_"+ye+mo+da+ho+mi+"_"+str(i).zfill(5)+'.png')
    plt.close()

#==========================================================
# End code
#==========================================================

