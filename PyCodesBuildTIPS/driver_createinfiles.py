def driver_createinfiles(i):
  "i corresponds to the imerg file"

  # Import python libraries
  import numpy as np
  from netCDF4 import Dataset
  import scipy.ndimage as ndimage
  import h5py
  import time

#==========================================================
# Read files output by main script
#==========================================================

  # Read in namelist
  nl = Dataset("namelist_ci.nc","r")

  # Read in filename
  fn = open("filenames_ci.txt")
  filenames = fn.read().split(',')

  print("Reading data from "+filenames[i])

  if nl.datanc4=="True":
    # Open file
    dataset = Dataset(filenames[i])

    # Read variables
    rain = dataset.variables["precipitationCal"][:,:,:]

  if nl.datahdf5=="True":
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
  lenIdir  = len(nl.datadirin)
  ye = str(filenames[i][lenIdir+21:lenIdir+25])
  mo = str(filenames[i][lenIdir+25:lenIdir+27])
  da = str(filenames[i][lenIdir+27:lenIdir+29])
  ho = str(filenames[i][lenIdir+31:lenIdir+33])
  mi = str(filenames[i][lenIdir+33:lenIdir+35])
  timenow = ye+mo+da+ho+mi

#==========================================================
# Subset data
#==========================================================

  if nl.ssreg=="True":

    print("Subsetting variable by area")

    if nl.datanc4=="True":
      # Read in coordinate variables
      lat = dataset.variables["lat"][:]
      lon = dataset.variables["lon"][:]

    if nl.datahdf5=="True":
      # Read in coordinate variables
      lat = datasetIM["Grid/lat"][:]
      lon = datasetIM["Grid/lon"][:]

    # Find indices of closest coordinates
    idyS = (np.abs(lat - nl.latS)).argmin()
    idyN = (np.abs(lat - nl.latN)).argmin()
    idxW = (np.abs(lon - nl.lonW)).argmin()
    idxE = (np.abs(lon - nl.lonE)).argmin()

    # Subset main dataset by indices
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

  if nl.smooth=="True":
    if nl.gaussian=="True":
      rain = ndimage.gaussian_filter(rain,sigma=nl.stddev)
    elif nl.uniform=="True":
      rain = ndimage.uniform_filter(rain,size=nl.width)
    
#==========================================================
# Do thresholding
#==========================================================

  print("Creating thresholded data")
  
  for j in (range(len(nl.tholds)+1)):
    if j==0:
      value1 = np.where(rain<nl.tholds[0],j,rain)
    elif j>0 and j<len(nl.tholds):
      value1 = np.where((rain>=nl.tholds[j-1]) & (rain<nl.tholds[j]),j,value1)
    elif j==len(nl.tholds):
      value1 = np.where(rain>=nl.tholds[j-1],j,value1)

#==========================================================
# Write thresholds to netcdf file
#==========================================================

  # Creating netcdf file to write to
  if nl.ssreg=="True": 
    print("Creating file: "+nl.datadirout+nl.fileidout+nl.ssname+"_"+str(i).zfill(5)+".nc")
    fileout = Dataset(nl.datadirout+nl.fileidout+nl.ssname+"_"+str(i).zfill(5)+".nc",'w',format='NETCDF4_CLASSIC')
  else:
    print("Creating file: "+nl.datadirout+nl.fileidout+str(i).zfill(5)+".nc")
    fileout = Dataset(nl.datadirout+nl.fileidout+str(i).zfill(5)+".nc",'w',format='NETCDF4_CLASSIC')

  # Global Attributes
  fileout.description = 'IMERG thresholds for FiT Tracking'
  fileout.history     = 'Created ' + time.ctime(time.time())
  fileout.source      = 'https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGHH_V06/summary?keywords=imerg'
  fileout.time        = str(timenow)
  fileout.tholds      = nl.tholds
  fileout.datestart   = str(nl.starttime)
  fileout.dateend     = str(nl.endtime)
  if nl.ssreg=="True":
    fileout.latN      = nl.latN
    fileout.latS      = nl.latS
    fileout.lonW      = nl.lonW
    fileout.lonE      = nl.lonE
    fileout.region    = nl.ssname
  else:
    fileout.region    ='global'
  if nl.smooth=="True":
    if nl.gaussian=="True":
      fileout.smthstddev=str(nl.stddev)  
    if nl.uniform=="True":
      fileout.smthwidth=str(nl.width)

  # Create dimensions in file
  z = fileout.createDimension('z', 1)
  y = fileout.createDimension('y', len(y))
  x = fileout.createDimension('x', len(x))

  # Create variables in file
  value  = fileout.createVariable('value' , np.float64, ('z','y','x'))

  # Write Variables
  value[:] = value1

#==========================================================
# End code
#==========================================================
