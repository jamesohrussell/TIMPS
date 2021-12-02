#==========================================================
# Script to read IMERG data files and convert them to files
# which can be read by the FiT tracking algorithm. Converts
# rain rates to thresholded data i.e. everything below a 
# low threshold is 0, everything between the first two 
# thresholds is 1, between second two thresholds is 2, 
# and above the last threshold is 3.
#
# Addition of Mani Rajagopal's normalized thresholds 
# methods (6/22/2019). These look at each contiguous 
# object and normalize the thresholds within each relative
# to the maximum precipation. Thus thresholds are a 
# fraction of the max precip within a contiguous object.
#
# James Russell 2019
#==========================================================

#==========================================================
# Import libraries
#==========================================================

# Import python libraries
import glob
import time
import datetime
from joblib import Parallel, delayed
import os

import sys
import numpy as np
from netCDF4 import Dataset
import scipy.ndimage as ndimage
import h5py

# Import namelist
from namelist_TIMPS import general as gnl
from namelist_TIMPS import create as cnl

# Import scripts
sys.path.insert(0,gnl.fnsdir)
import shape_functions as sfns
import misc_functions as mfns
import time_functions as tfns

#==========================================================
# Initialize timer
#==========================================================

startnow = time.time()

#==========================================================
# Read files
#==========================================================

print("Generating file list")

# Generate a list of filenames with dates to search for
start = datetime.datetime.strptime(gnl.starttime,"%Y%m%d")
end = datetime.datetime.strptime(gnl.endtime,"%Y%m%d")
datearr = (start+datetime.timedelta(days=x) 
 for x in range(0,(end-start).days))
filen = []

#!! Change this line for different directory structure !!#
for dateobj in datearr:
  filen.append(f'{gnl.datadirin}{dateobj.strftime("%Y")}/\
{dateobj.strftime("%m")}/{gnl.fileidin}\
{dateobj.strftime("%Y%m%d")}*')

# Reads directory and filenames
n = 0
for f in filen:
  filestoadd = glob.glob(f)

  # All half hourly files must be present to continue
  # Ensures the time-series is continuous
  if not len(filestoadd)==48:
    if len(filestoadd)<48:
      raise ValueError(f"Missing files for {f}")
    else:
      raise ValueError(f"Too many files for {f}")

  # Adds new files to list
  if n==0:
    filenames = filestoadd
  else:
    filenames = filenames+filestoadd
  n=n+1

if len(filenames)==0:
  raise ValueError("No IMERG files found. Perhaps the IMERG \
directory or fileid is incorrect. Please check the \
namelist and try again.")

#==========================================================
# Define function to run
#==========================================================

def driver_createinfiles(i):
  "i corresponds to the file"

#==========================================================
# Read files output by main script

  print(f'Reading data from {filenames[i]}')

  if gnl.datanc4:
    # Open file
    dataset = Dataset(filenames[i])

    # Read variables
    rain = dataset.variables["precipitationCal"][:,:,:]
    tmIM = dataset.variables["time"][:]
    tmun = dataset.variables["time"].getncattr("units")

  if gnl.datahdf5:
    # Open file
    datasetIM = h5py.File(filenames[i],'r')
    # Read variables
    rain = datasetIM["/Grid/precipitationCal"][:,:,:]
    tmIM = datasetIM["/Grid/time"][:]
    tmun = datasetIM["/Grid/time"].attrs.get("units")

  # Change dimension order of variable
  rain = np.transpose(rain, (0, 2, 1))

  # Get start time in YYYYMMDDhhmm format
  tn1 = tfns.time_since_inv(
   float(tmIM[0]),tmun.decode("UTF-8"))
  timenow = tn1[0:4]+tn1[5:7]+tn1[8:10]+\
   tn1[11:13]+tn1[14:16]

#==========================================================
# Subset data

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
    lat  = lat[idyS:idyN+1]
    lon  = lon[idxW:idxE+1]
    rain = rain[:,idyS:idyN+1,idxW:idxE+1]

  # Get subset dimensions and set x and y coordinates
  szrain = np.shape(rain)
  x = np.linspace(0, szrain[2]-1, szrain[2])
  y = np.linspace(0, szrain[1]-1, szrain[1])

#==========================================================
# Smooth data

  if cnl.smooth:
    if cnl.gaussian:
      rain = ndimage.gaussian_filter(rain,sigma=cnl.stddev)
    elif cnl.uniform:
      rain = ndimage.uniform_filter(rain,size=cnl.width)
    
#==========================================================
# Do thresholding

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
    #objmax = [10 if i<10 else i for i in objmax]
    rainobjmax = np.copy(rainlabs)
    for il in labels:
      rainobjmax = np.where(rainobjmax==il,
         objmax[il-1],rainobjmax)
    
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

  # Contiguous area
  elif cnl.tholdtype==3:

    # Set everything above min threhold to zero
    value1 = np.where(rain>=cnl.minthold,1,0)

  else:
    raise ValueError("Incorrect option for threshold")

#==========================================================
# Write thresholds to netcdf file

  # Creating netcdf file to write to
  if cnl.ssreg: 
    print(f'Creating file: {gnl.datadirFiTin}\
{gnl.fileidFiTin}{str(i).zfill(5)}.nc')
    fileout = Dataset(f'{gnl.datadirFiTin}\
{gnl.fileidFiTin}{str(i).zfill(5)}.nc','w',
     format='NETCDF4_CLASSIC')
  else:
    print(f'Creating file: {gnl.datadirFiTin}\
{gnl.fileidFiTin}{str(i).zfill(5)}.nc')
    fileout = Dataset(gnl.datadirFiTin+gnl.fileidFiTin+ \
     str(i).zfill(5)+".nc",'w',format='NETCDF4_CLASSIC')

  # Global Attributes
  fileout.description = 'IMERG thresholds for FiT Tracking'
  fileout.history     = f'Created {time.ctime(time.time())}'
  fileout.source      = 'https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGHH_V06/summary?keywords=imerg'
  fileout.time        = str(timenow)
  fileout.tholds      = cnl.tholds
  fileout.datestart   = str(gnl.starttime)
  fileout.dateend     = str(gnl.endtime)
  if cnl.ssreg:
    fileout.latN      = max(lat)
    fileout.latS      = min(lat)
    fileout.lonW      = min(lon)
    fileout.lonE      = max(lon)
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
# End function
#==========================================================

#==========================================================
# Parallel loop over object
#========================================================== 

# Begins loop
if cnl.serialorparallel==1:
  print("Begin serial loop over objects")
  for i in range(len(filenames)):
    driver_createinfiles(i)

# Parallel loop over PFs
if cnl.serialorparallel==2:
  print("Begin parallel loop over objects")
  Parallel(n_jobs=cnl.njobs)(delayed(
    driver_createinfiles)(i) for \
    i in range(len(filenames)))

#==================================================================
# Generate namelist to run tracking
#==================================================================

ffn = open(f'{gnl.datadirFiTin}namelist_FiT.txt',"w")
ffn.write("input_file_type = 2\n")
ffn.write("output_file_type = 2\n")
ffn.write(f'input_file_prefix = {gnl.datadirFiTin}\
{gnl.fileidFiTin}\n')
ffn.write(f'output_file_prefix = {gnl.datadirFiTout}\
{gnl.fileidFiTout1}\n')
file1 = Dataset(f'{gnl.datadirFiTin}{gnl.fileidFiTin}00000.nc')
ffn.write(f'dimx = {str(file1.dimensions["x"].size)}\n')
ffn.write(f'dimy = {str(file1.dimensions["y"].size)}\n')
file1.close()
ffn.write("dimz = 1\n")
ffn.write(f'number_of_timesteps = {str(len(filenames))}\n')
if cnl.tholdtype==3:
  ffn.write("number_of_thresholds = 1\n")
else:
  ffn.write(f'number_of_thresholds = \
{str(len(cnl.tholds)+1)}\n')
ffn.write(f'horizontal_distance_limit_for_merging = \
{str(gnl.splitdist)}\n')
if gnl.domperiodicx:
  ffn.write("domain_periodic_in_x_dimension = 1\n")
else:
  ffn.write("domain_periodic_in_x_dimension = 0\n")
ffn.write("output_3D_objects = 0;\n")
ffn.write("output_4D_objects = 1;\n")
ffn.write("output_png_files = 0;\n")

#==================================================================
# Final clean up tasks
#==================================================================

# Remove namelist and filename files
print("Cleaning up")
os.system("rm -rf __pycache__")

# End timer
endnow = time.time()
print(f'Program took {str(endnow - startnow)}s')

#==================================================================
# End code
#==================================================================

