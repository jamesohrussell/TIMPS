# Import python libraries
import numpy as np
from netCDF4 import Dataset
import scipy.ndimage as ndimage
import h5py
import time
import sys

# Import namelist
from namelist_TIPS import general as gnl
from namelist_TIPS import create as cnl

# Import custom libraries
sys.path.insert(0,gnl.fnsdir)
import shape_functions as sfns
import misc_functions as mfns

def driver_createinfiles(i):
  "i corresponds to the file"

#==========================================================
# Read files output by main script
#==========================================================

  # Read in filename
  fn = open("filenames_ci.txt")
  filenames = fn.read().split(',')

  print("Reading data from "+filenames[i])

  if gnl.datanc4:
    # Open file
    dataset = Dataset(filenames[i])

    # Read variables
    rain = dataset.variables["precipitationCal"][:,:,:]

  if gnl.datahdf5:
    # Open file
    datasetIM = h5py.File(filenames[i],'r')

    # Read variables
    rain = datasetIM["/Grid/precipitationCal"][:,:,:]

  # Change dimension order of variable
  rain = np.transpose(rain, (0, 2, 1))

  # Obtain time and date, compile into string in format:
  #   YYYYMMDDhhmm
  lenIdir  = len(gnl.datadirin)
  ye = str(filenames[i][lenIdir+21:lenIdir+25])
  mo = str(filenames[i][lenIdir+25:lenIdir+27])
  da = str(filenames[i][lenIdir+27:lenIdir+29])
  ho = str(filenames[i][lenIdir+31:lenIdir+33])
  mi = str(filenames[i][lenIdir+33:lenIdir+35])
  timenow = ye+mo+da+ho+mi

#==========================================================
# Subset data
#==========================================================

  if cnl.ssreg:

    print("Subsetting variable by area")

    if gnl.datanc4:
      # Read in coordinate variables
      lat = dataset.variables["lat"][:]
      lon = dataset.variables["lon"][:]

    if gnl.datahdf5:
      # Read in coordinate variables
      lat = datasetIM["Grid/lat"][:]
      lon = datasetIM["Grid/lon"][:]

    # Find indices of closest coordinates
    idyS = (np.abs(lat - cnl.latS)).argmin()
    idyN = (np.abs(lat - cnl.latN)).argmin()
    idxW = (np.abs(lon - cnl.lonW)).argmin()
    idxE = (np.abs(lon - cnl.lonE)).argmin()

    # Subset main dataset by indices
    lat  = lat[idyS:idyN]
    lon  = lon[idxW:idxE]
    rain = rain[:,idyS:idyN,idxW:idxE]

  # Get subset dimensions and set x and y coordinates
  szrain = np.shape(rain)
  x = np.linspace(0, szrain[2]-1, szrain[2])
  y = np.linspace(0, szrain[1]-1, szrain[1])

#==========================================================
# Smooth data
#==========================================================

  if cnl.smooth:
    if cnl.gaussian:
      rain = ndimage.gaussian_filter(rain,sigma=cnl.stddev)
    elif cnl.uniform:
      rain = ndimage.uniform_filter(rain,size=cnl.width)
    
#==========================================================
# Do thresholding
#==========================================================

  print("Creating thresholded data")

  # Constant thresholds
  if cnl.tholdtype==1:
  
    # Loop over all thresholds
    for j in (range(len(cnl.tholds)+1)):
      if j==0:
        value1 = np.where(rain<cnl.tholds[0],j,rain)
      elif j>0 and j<len(cnl.tholds):
        value1 = np.where((rain>=cnl.tholds[j-1]) & \
                          (rain<cnl.tholds[j]),j,value1)
      elif j==len(cnl.tholds):
        value1 = np.where(rain>=cnl.tholds[j-1],j,value1)

  # Normalized thresholds
  elif cnl.tholdtype==2:
    
    # Remove all precip below minthold
    rain = np.where(rain>=cnl.minthold,rain-cnl.minthold,0)

    # Find all contiguous objects
    rainlabs,numobs = sfns.label_wdiags(rain)
    labels = np.arange(1,numobs+1)

    # Get max information for each object and
    # make new array with max of each object
    objmax = ndimage.labeled_comprehension(
     rain,rainlabs,labels,np.max,np.double,0)
    rainobjmax = np.copy(rainlabs)
    for il in labels:
      rainobjmax[rainlabs==il] = objmax[il-1]
    
    # Normalize precip
    nrain = np.divide(rain,rainobjmax,
     out=np.zeros_like(rain),where=rainobjmax!=0)    

    # Do thresholding
    for j in (range(len(cnl.tholds)+1)):
      if j==0:
        value1 = np.where((nrain>0.) & \
         (nrain<cnl.tholds[0]),j+1,nrain)
      elif j>0 and j<len(cnl.tholds):
        value1 = np.where((nrain>=cnl.tholds[j-1]) & \
                          (nrain<cnl.tholds[j]),j+1,value1)
      elif j==len(cnl.tholds):
        value1 = np.where(nrain>=cnl.tholds[j-1],j+1,value1)

  else:
    raise ValueError("Incorrect option for threshold")

#==========================================================
# Write thresholds to netcdf file
#==========================================================

  # Creating netcdf file to write to
  if cnl.ssreg: 
    print("Creating file: "+cnl.datadirout+cnl.fileidout+ \
     cnl.ssname+"_"+str(i).zfill(5)+".nc")
    fileout = Dataset(cnl.datadirout+cnl.fileidout+ \
     cnl.ssname+"_"+str(i).zfill(5)+".nc",'w',
     format='NETCDF4_CLASSIC')
  else:
    print("Creating file: "+cnl.datadirout+ \
     cnl.fileidout+str(i).zfill(5)+".nc")
    fileout = Dataset(cnl.datadirout+cnl.fileidout+ \
     str(i).zfill(5)+".nc",'w',format='NETCDF4_CLASSIC')

  # Global Attributes
  fileout.description = 'IMERG thresholds for FiT Tracking'
  fileout.history     = 'Created ' + time.ctime(time.time())
  fileout.source      = 'https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGHH_V06/summary?keywords=imerg'
  fileout.time        = str(timenow)
  fileout.tholds      = cnl.tholds
  fileout.datestart   = str(cnl.starttime)
  fileout.dateend     = str(cnl.endtime)
  if cnl.ssreg:
    fileout.latN      = cnl.latN
    fileout.latS      = cnl.latS
    fileout.lonW      = cnl.lonW
    fileout.lonE      = cnl.lonE
    fileout.region    = cnl.ssname
  else:
    fileout.region    ='global'
  if cnl.smooth:
    if cnl.gaussian:
      fileout.smthstddev=str(cnl.stddev)  
    if cnl.uniform:
      fileout.smthwidth=str(cnl.width)

  # Create dimensions in file
  z = fileout.createDimension('z', 1)
  y = fileout.createDimension('y', len(y))
  x = fileout.createDimension('x', len(x))

  # Create variables in file
  value  = fileout.createVariable('value',
           np.float64, ('z','y','x'))

  # Write Variables
  value[:] = value1

#==========================================================
# End code
#==========================================================
