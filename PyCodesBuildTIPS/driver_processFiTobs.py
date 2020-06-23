#============================================================
# Import libraries
#============================================================

# Import python libraries
import numpy as np
import h5py
import csv
import scipy.ndimage as ndimage
import re
from netCDF4 import Dataset
import datetime as dt
import time as tm
import sys

# Import namelist
from namelist_TIPS import general as gnl
from namelist_TIPS import process as pnl

# Load custom libraries
sys.path.insert(0,gnl.fnsdir)
import misc_functions as mfns

#============================================================
# Begin function
#============================================================

def driver_processFiTobs(o):
  "o corresponds to the object in objs"

#============================================================
# Read files output by main script
#============================================================

  # Read text file
  f       = open(pnl.datadirFi+pnl.fileidFi1+pnl.fileidtxt,'r') 
  FiTinfo = np.genfromtxt(f, dtype="str",delimiter="\t",
                           skip_header=1)

  # Open netcdf file and assign data
  objs = Dataset("passinfo_po.nc").objs

  # Read in filename
  fn = open("filenames_po.txt")
  filenames = fn.read().split(',')

#============================================================
# Begin loop over times for object
#============================================================
  
  # Find all time indices and times relevant for object
  indob = [i for i, x in enumerate(list(
             map(int,FiTinfo[:,0]))) if x == int(objs[o])]
  times = list(map(int,FiTinfo[indob,2]))

  print("Working on object "+str(objs[o])+" with "+ \
         str(len(times))+" times")

  # Times in hours since reftime
  time1          = [-999.]* len(times) 
  # Times (YYYYMMDDhhmm)
  datetime       = [-999] * len(times)
  # Central latitude of object
  centrallat     = [-999.] * len(times) 
  # Central longitude of object
  centrallon     = [-999.] * len(times)
  # Weighted central latitude of object
  wgtcentlat     = [-999.] * len(times) 
  # Weighted central longitude of object
  wgtcentlon     = [-999.] * len(times)
  # Dictionaries
  lats = {}; lons = {}; instrain = {}; QI = {}

  # Set reference time
  sincedate = dt.datetime.strptime(pnl.reftime, "%Y-%m-%d %H:%M:%S")

  # Loop over all times relevant to object
  for t in range(0,len(times)):

#============================================================
# Read FiT objects file and process for specific object
#============================================================

    # Open FiT input file
    Finfilename = pnl.datadirFin+pnl.fileidFin+str(
                    times[t]).zfill(5)+".nc"
    datasetFin  = Dataset(Finfilename)
    timenow     = str(datasetFin.time)

    # Obtain universal information
    region = datasetFin.region
    tholds = datasetFin.tholds
    source = datasetFin.source
    if region=='global':
      ssreg    = False
    else:
      ssreg    = True
      latN     = datasetFin.latN
      latS     = datasetFin.latS
      lonW     = datasetFin.lonW
      lonE     = datasetFin.lonE
    if hasattr(datasetFin, 'smthstddev'):
      sstd     = datasetFin.smthstddev
    else:
      sstd     = None
    if hasattr(datasetFin, 'smthwidth'):
      swid     = datasetFin.smthwidth
    else:
      swid     = None

    # Open file and read map of object ids
    FiTfilename = pnl.datadirFi+pnl.fileidFi1+str(
                    times[t]).zfill(5)+pnl.fileidFi2
    datasetFi   = Dataset(FiTfilename)
    objects     = np.squeeze(
                    datasetFi.variables["value"][:,:,:])

    # Make obbin only 1 for specific object id
    obbin       = np.where(objects==int(objs[o]),1,0)

    # Find all location indices of object
    locindob    = np.where(obbin==1)

#============================================================
# Read IMERG data and process for specific object
#============================================================
   
    if gnl.datanc4:
      # Open file
      datasetIM = Dataset(filenames[times[t]])

      # Read variables
      rain  = datasetIM.variables["precipitationCal"][:,:,:]
      QI1  = datasetIM.variables[
       "precipitationQualityIndex"][:,:,:]

    if gnl.datahdf5:
      # Open file
      datasetIM = h5py.File(filenames[times[t]],'r')
      # Read variables
      rain = datasetIM["/Grid/precipitationCal"][:,:,:]
      QI1 = datasetIM[
       "Grid/precipitationQualityIndex"][:,:,:]

    # Change dimension order of variable
    rain = np.squeeze(np.transpose(rain,(0, 2, 1)))
    QI1  = np.squeeze(np.transpose(QI1,(0, 2, 1)))

    # Subset data
    if ssreg:

      # Read in coordinate variables
      if gnl.datanc4:
        lat = datasetIM.variables["lat"][:]
        lon = datasetIM.variables["lon"][:]
      
      if gnl.datahdf5:
        lat = datasetIM["Grid/lat"][:]
        lon = datasetIM["Grid/lon"][:]
        
      # Find indices of closest coordinates
      idyS = (np.abs(lat - latS)).argmin()
      idyN = (np.abs(lat - latN)).argmin()
      idxW = (np.abs(lon - lonW)).argmin()
      idxE = (np.abs(lon - lonE)).argmin()

      # Subset main dataset by indices
      lat  = lat[idyS:idyN]
      lon  = lon[idxW:idxE]
      rain = rain[idyS:idyN,idxW:idxE]
      QI1  = QI1[idyS:idyN,idxW:idxE]

#============================================================
# Calculate parameters and assign data at first time
#============================================================
    
    # Assign some simple values
    datetime[t]    = int(timenow)
    centrallat[t]  = np.mean(lat[locindob[0]])
    centrallon[t]  = np.mean(lon[locindob[1]])

    # Calculate weighted (by rainfall) centroid
    sr = sum(rain[locindob[0],locindob[1]])
    if sr!=0: 
      rw = [mfns.divzero(r,sr) for r in rain[locindob[0],locindob[1]]]
      wgtcentlat[t] = np.average(lat[locindob[0]],weights=rw)
      wgtcentlon[t] = np.average(lon[locindob[1]],weights=rw)

    # Make reference time object and calculate time
    #  in hours since reftime
    currentdate = dt.datetime.strptime(timenow,"%Y%m%d%H%M")
    delta = (currentdate - sincedate)
    time1[t] = delta.total_seconds()/3600.

    # For each new time, add to dictionary as a new key
    lats[str(datetime[t])] = lat[locindob[0]]
    lons[str(datetime[t])] = lon[locindob[1]]
    instrain[str(datetime[t])] = \
     rain[locindob[0],locindob[1]]
    QI[str(datetime[t])] = QI1[locindob[0],locindob[1]]

  # For PF, if all rain rates at all times are actually zero
  #  move onto the next PF
  for i in instrain.keys():
    if len(instrain[i][instrain[i]>0])==0:
      a = True
    else:
      a = False
      break
  if a:
    return
  
  # For all other PFs, if all rain rates at any given time 
  #  are zero, add that index to a list of indices to delete
  delind = []
  ind = 0
  for i in instrain.keys():
    if len(instrain[i][instrain[i]>0])==0:
      delind.append(ind)
    ind = ind+1

  # Delete the key, value pairs in the dictionaries
  for i in delind:
    del(lats[str(datetime[i])])
    del(lons[str(datetime[i])]) 
    del(instrain[str(datetime[i])])
    del(QI[str(datetime[i])])

  # Delete the indices in the lists
  delind.reverse()
  for i in delind:
    del(datetime[i])
    del(time1[i])
    del(centrallat[i])
    del(centrallon[i])
    del(wgtcentlat[i])
    del(wgtcentlon[i])

#============================================================
# Write NetCDF file for object
#============================================================

  # Creating netcdf file to write to
  filename = gnl.fileidTIPS+str(objs[o])+"_"+ \
   str(datetime[0])+"_"+ \
   str(int(round(centrallat[0]))).zfill(2)+"_"+ \
   str(int(round(centrallon[0]))).zfill(2)+".nc"

  print("Creating file: "+gnl.datadirTIPS+filename)
  fileout = Dataset(gnl.datadirTIPS+filename,
   'w',format='NETCDF4')

  # Global Attributes
  fileout.long_name   = 'Precipitation feature information'
  fileout.description = 'Information for a precipitation feature tracked by FiT algorithm with specifications as in the attributes below'
  fileout.history     = 'Created ' + tm.ctime(tm.time())
  fileout.source      = 'IMERG V06B: '+source
  fileout.objectID    = str(objs[o])
  if ssreg:
    fileout.track_dom_latN   = latN
    fileout.track_dom_latS   = latS
    fileout.track_dom_lonW   = lonW
    fileout.track_dom_lonE   = lonE
    fileout.track_dom_region = region 
  if sstd is not None:
    fileout.smthcommand = "ndimage.gaussian_filter(rain,sigma=int(smthstdd))"
    fileout.smthstdd = stdd
  if swid is not None:
    fileout.smthcommand = "ndimage.uniform_filter(rain,size=int(smthwidth))"
    fileout.smthwidth = swid
    fileout.smthwidth_units  = "pixels"
    fileout.smth_description = "Smoothing used for tracking objects in FiT algorithm"
  fileout.splitdist   = int(pnl.mergedist)
  fileout.splitdist_units = "pixels"
  fileout.splitdist_description = "Maximum distance objects can be separated but remain the same object after splitting in the FiT algorithm"
  fileout.tholds      = tholds
  fileout.tholds_units = "mm/hr"
  fileout.tholds_description = "Rain rate thresholds used for tracking with FiT algorithm"

  # Create dimensions in file
  time = fileout.createDimension('time', len(time1))

  # Create variables in file
  mfns.write_var("time","Time","","time",np.float64,
                 "hours since "+pnl.reftime,fileout,
                 time1,gnl.datadirTIPS+filename,-999.)

  mfns.write_var("datetime","Date and time","","time",
                 np.int64,"YYYYMMDDhhmm",fileout,datetime,
                 gnl.datadirTIPS+filename,-999)

  mfns.write_var("centrallat","Central latitude",
                 "Latitude of PF centroid","time",
                 np.float64,"degreesNorth",fileout,
                 centrallat,gnl.datadirTIPS+filename,-999.)

  mfns.write_var("centrallon","Central longitude",
                 "Longitude of PF centroid","time",
                 np.float64,"degreesEast",fileout,
                 centrallon,gnl.datadirTIPS+filename,-999.)

  description = "Latitude of PF centroid weighted by rainfall"
  mfns.write_var("centlatwgt","Weighted Central Latitude",
                 description,"time",np.float64,"degreesNorth",fileout,
                 wgtcentlat,gnl.datadirTIPS+filename,-999.)

  description = "Longitude of PF centroid weighted by rainfall"
  mfns.write_var("centlonwgt","Weighted Central Longitude",
                 description,"time",np.float64,"degreesEast",fileout,
                 wgtcentlon,gnl.datadirTIPS+filename,-999.)

  format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time."

  mfns.write_group("lats","Latitudes",
                 "Latitudes of IMERG grid cell centers for PF",
                 "DegreesNorth",format1,fileout,
                 lats,gnl.datadirTIPS+filename)

  mfns.write_group("lons","Longitudes",
                 "Longitudes of IMERG grid cell centers for PF",
                 "DegreesEast",format1,fileout,
                 lons,gnl.datadirTIPS+filename)

  description = "Instantaneous rain rates at IMERG grid cells corresponding to all latitude and longitude values in PF."
  mfns.write_group("instrain","Instantaneous rain rate",
                 description,"mm/hr",format1,fileout,
                 instrain,gnl.datadirTIPS+filename)

  description = "Quality index for rain rates at IMERG grid cells corresponding to all latitude and longitude values in PF. 0 is bad. 1 is good. Rule of thumb from IMERG team: Below 0.3 is bad."
  mfns.write_group("QI","Quality index",
                 description,"mm/hr",format1,fileout,
                 QI,gnl.datadirTIPS+filename)

  # Close file
  fileout.close()

