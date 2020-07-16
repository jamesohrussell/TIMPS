# Import libraries
import scipy.ndimage as ndimage
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps
import h5py
from collections import Counter
import csv
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

# Import namelist
import namelist_plot as nl

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
# Function to read text file
#==========================================================

def driver_plotIMERG(i):
  "Driver for plotting IMERG and FiT objects/thresholds in parrallel"

  # Read in filename
  for fn, row in enumerate(open("filenames_plot.txt")):
    if fn==i:
      f = row[:-1]

  print("Working on "+f)

  lenIdir = len(nl.datadirIM)

#==========================================================
# Read IMERG data (NetCDF4)
#==========================================================

  if nl.datanc4:

    # Open file
    datasetIM = Dataset(f)

    # Read variables
    rain = datasetIM.variables["precipitationCal"][:,:,:]

#==========================================================
# Read IMERG data (hdf5)
#==========================================================

  if nl.datahdf5:

    # Open file
    print(f)
    datasetIM = h5py.File(f,'r')

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

  if nl.plotobj:
    # Open file
    datasetFi = Dataset(nl.datadirFi+nl.fileidFi1+str(i).zfill(5)+nl.fileidFi2)
  
    # Read variables
    obid = np.squeeze(datasetFi.variables["value"][:,:,:])

#==========================================================
# Process FiT tracking data
#==========================================================

    if nl.idobj:
      # Make obbin only 1 for specific object id
      obbin = np.where(obid==obid1,1,0)
    elif nl.longobj:
      # Select largest nlngst objects
      datatxt = read_FiT_txt(nl.datadirFi+nl.fileidFi3)
      obidsdrs = Counter([int(i) for i in datatxt[:,0]])
      obidssrt = sorted(obidsdrs, key=obidsdrs.get,reverse=True)
      obids = obidssrt[0:nlong]

#==========================================================
# Take time from IMERG file name
#==========================================================
  
  ye = str(f[lenIdir+21:lenIdir+25])
  mo = str(f[lenIdir+25:lenIdir+27])
  da = str(f[lenIdir+27:lenIdir+29])
  ho = str(f[lenIdir+31:lenIdir+33])
  mi = str(f[lenIdir+33:lenIdir+35])

#==========================================================
# Subset data
#==========================================================

  if nl.ssreg:

    print("Subsetting variable by area")
    
    # Read in coordinate variables
    if nl.datanc4:
      lat = datasetIM.variables["lat"][:]
      lon = datasetIM.variables["lon"][:]

    if nl.datahdf5:
      lat = datasetIM["Grid/lat"][:]
      lon = datasetIM["Grid/lon"][:]

    # Find indices of closest coordinates
    idyS = (np.abs(lat - nl.latS)).argmin()
    idyN = (np.abs(lat - nl.latN)).argmin()
    idxW = (np.abs(lon - nl.lonW)).argmin()
    idxE = (np.abs(lon - nl.lonE)).argmin()

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

  if nl.smooth:
    print("Smoothing rainfall")
    if nl.gaussian:
      rain = ndimage.gaussian_filter(rain, sigma=nl.stddev)
    elif nl.uniform:
      rain = ndimage.uniform_filter(rain,size=nl.width)

#==========================================================
# Make plot
#==========================================================

  if nl.plot:
    print("Plotting data")

    # Generate center points for contourf
    lonc,latc = np.meshgrid(lon,lat) 

    # Generate point bounds for pcolorf
    lats, lons = np.mgrid[slice(nl.latS, nl.latN + nl.dy, nl.dy),
                    slice(nl.lonW, nl.lonE + nl.dx, nl.dx)]

    # Make figure
    ax = plt.figure(num=None, figsize=(8, 6), dpi=200, facecolor='w', edgecolor='k')
    #ax = plt.subplot()
  
    # set up a map
    extent = [nl.lonW, nl.lonE, nl.latS, nl.latN]
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(extent,crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m')

    # Format tick marks
    ax.set_xticks(np.linspace(nl.lonW,nl.lonE,num=nl.labnx), crs=ccrs.PlateCarree())
    ax.set_yticks(np.linspace(nl.latS,nl.latN,num=nl.labny), crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    # Make plot
    cs = plt.pcolormesh(lons,lats,np.squeeze(rain),\
                        transform=ccrs.PlateCarree(),\
                        cmap=cmaps.prcp_1,\
                        vmin=nl.mincl,vmax=nl.maxcl)
    cbar = plt.colorbar(cs, ax=ax, orientation="horizontal", pad=.1, fraction=0.06)   
    cbar.set_label("IMERG Precipitation (mm/hr)")

    # Contour outline of objects
    if nl.plotobj:

      cntrclr = ['r','g','b','m','y','#1f77b4','#ff7f0e','#9467bd','#8c564b','#7f7f7f']        
     
      if nl.longobj:
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
    if nl.plotthr:
      plt.contour(lonc, latc, np.squeeze(rain), nl.tholds, colors='k')

    plt.xlim(nl.lonmn, nl.lonmx)
    plt.ylim(nl.latmn, nl.latmx)
    plt.annotate(ye+"/"+mo+"/"+da+" "+ho+":"+mi, xy=(0.78, 0.92), xycoords='axes fraction')
    plt.savefig('figs/IMERG_'+nl.ssname+"_"+ye+mo+da+ho+mi+"_"+str(i).zfill(5)+'.png')
    plt.close()

#==========================================================
# End code
#==========================================================

