#==========================================================
# Script to read IMERG data files and convert them to files
# which can be read by the FiT tracking algorithm. Converts
# rain rates to sort-of binary i.e. everything below a 
# low threshold is 0, everything between the first two 
# thresholds is 1, between second two thresholds is 2, 
# and above the last threshold is 3.
#
# James Russell 2019
#==========================================================

# Import python libraries
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import glob
import time
import scipy.ndimage as ndimage
import h5py
import datetime

#==========================================================
# Namelist
#==========================================================

# Directory and filename for input IMERG data
datadirin = "/uufs/chpc.utah.edu/common/home/varble-group2/IMERG/"
fileidin  = "3B-HHR.MS.MRG.3IMERG."
datahdf5  = True
datanc4   = False

# Date and time range
starttime = "20180601"
endtime   = "20180610" # Actual last day is day before

# Thresholds for rain rates (lowest to highest)
tholds  = [0.8,3.,9.,27.]

# Directory and filename for output FiT input data
datadirout = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_102519/FiT_input_files_2018_v2/"
fileidout  = "IMERG_FiT_tholds_"

# Subset regions (ranges have no affect if ssreg=False)
ssreg  = True
latN   = 30
latS   = -10
lonW   = -70
lonE   = 0
ssname = "CapVer"

# Make plot (for debugging)
plot = False

# Smoothing
smooth   = True
gaussian = False
stddev   = 1.5  # Standard deviation for gaussian smoothing
uniform  = True
width    = 5    # Number of points to spread running average

#==========================================================
# Initialize timer
#==========================================================

startnow = time.time()

#==========================================================
# Read data
#==========================================================

print("Generating list of filenames")

# Generate a list of filenames with dates to search for
start = datetime.datetime.strptime(starttime, "%Y%m%d")
end = datetime.datetime.strptime(endtime, "%Y%m%d")
datearr = (start + datetime.timedelta(days=x) for x in range(0, (end-start).days))
filen = []
for dateobj in datearr:
  filen.append(datadirin+fileidin+dateobj.strftime("%Y%m%d")+"*")

# Reads directory and filenames
n = 0
for f in filen:
  if n==0:
    filenames = glob.glob(f)
  else:
    filenames = filenames+glob.glob(f)
  n=n+1

# Begins loop
lenIdir  = len(datadirin)
for i in range(len(filenames)):

  print("Reading data from "+filenames[i])

  if datanc4:
    # Open file
    dataset = Dataset(filenames[i])

    # Read variables
    rain = dataset.variables["precipitationCal"][:,:,:]

  if datahdf5:
    # Open file
    datasetIM = h5py.File(filenames[i],'r')

    # Read variables
    rain = datasetIM["/Grid/precipitationCal"][:,:,:]

  # Change dimension order of variable
  raint = np.transpose(rain, (0, 2, 1))
  del(rain)
  rain = raint
  del(raint)

  # Obtain time and date, compile into string in format:
  #   YYYYMMDDhhmm
  ye = str(filenames[i][lenIdir+21:lenIdir+25])
  mo = str(filenames[i][lenIdir+25:lenIdir+27])
  da = str(filenames[i][lenIdir+27:lenIdir+29])
  ho = str(filenames[i][lenIdir+31:lenIdir+33])
  mi = str(filenames[i][lenIdir+33:lenIdir+35])
  timenow = ye+mo+da+ho+mi

#==========================================================
# Subset data
#==========================================================

  if ssreg:

    print("Subsetting variable by area")

    if datanc4:
      # Read in coordinate variables
      lat = dataset.variables["lat"][:]
      lon = dataset.variables["lon"][:]

    if datahdf5:
      # Read in coordinate variables
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
    del(rain,lat,lon)
    rain = rainsub
    lat = latsub
    lon = lonsub
    del(rainsub,latsub,lonsub)

  # Get subset dimensions and set x and y coordinates
  szrain = np.shape(rain)
  x = np.linspace(0, szrain[2]-1, szrain[2])
  y = np.linspace(0, szrain[1]-1, szrain[1])

#==========================================================
# Smooth data
#==========================================================

  if smooth:
    if gaussian:
      rain = ndimage.gaussian_filter(rain,sigma=stddev)
    elif uniform:
      rain = ndimage.uniform_filter(rain,size=width)
    
#==========================================================
# Do thresholding
#==========================================================

  print("Creating thresholded data")
  
  for j in (range(len(tholds)+1)):
    if j==0:
      value1 = np.where(rain<tholds[0],j,rain)
    elif j>0 and j<len(tholds):
      value1 = np.where((rain>=tholds[j-1]) & (rain<tholds[j]),j,value1)
    elif j==len(tholds):
      value1 = np.where(rain>=tholds[j-1],j,value1)

#==========================================================
# Make plot to confirm data is correctly read
#==========================================================

  if plot:
    X,Y = np.meshgrid(x,y)
    plt.contourf(np.squeeze(X),np.squeeze(Y),np.squeeze(value1[0,:,:]))
    plt.show()

#==========================================================
# Write thresholds to netcdf file
#==========================================================

  # Creating netcdf file to write to
  if ssreg: 
    print("Creating file: "+datadirout+fileidout+ssname+"_"+str(i).zfill(5)+".nc")
    fileout = Dataset(datadirout+fileidout+ssname+"_"+str(i).zfill(5)+".nc",'w',format='NETCDF4_CLASSIC')
  else:
    print("Creating file: "+datadirout+fileidout+str(i).zfill(5)+".nc")
    fileout = Dataset(datadirout+fileidout+str(i).zfill(5)+".nc",'w',format='NETCDF4_CLASSIC')

  # Global Attributes
  fileout.description = 'IMERG thresholds for FiT Tracking'
  fileout.history     = 'Created ' + time.ctime(time.time())
  fileout.source      = 'https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGHH_V06/summary?keywords=imerg'
  fileout.time        = str(timenow)
  fileout.tholds      = tholds
  fileout.datestart   = str(starttime)
  fileout.dateend     = str(endtime)
  if ssreg:
    fileout.latN      = latN
    fileout.latS      = latS
    fileout.lonW      = lonW
    fileout.lonE      = lonE
    fileout.region    = ssname
  else:
    fileout.region    ='global'
  if smooth:
    if gaussian:
      fileout.smthstddev=str(stddev)  
    if uniform:
      fileout.smthwidth=str(width)

  # Create dimensions in file
  z = fileout.createDimension('z', 1)
  y = fileout.createDimension('y', len(y))
  x = fileout.createDimension('x', len(x))

  # Create variables in file
  value  = fileout.createVariable('value' , np.float64, ('z','y','x'))

  # Write Variables
  value[:] = value1


#==========================================================
# Final clean up tasks
#==========================================================

# End timer
endnow = time.time()
print("Program took "+str(endnow - startnow)+"s")

#==========================================================
# End code
#==========================================================
