#==========================================================
# Compile IMERG and FiT information into precipitation 
# feature netcdf files.
#
# Must provide in namelist:
#  * Directory and filename identifiers for IMERG files 
#      and FiT output netcdf files (the latter has no 
#      effect if only plotting IMERG)
#
# James Russell 2019
#==========================================================

# Import python libraries (do not change)
import numpy as np
from netCDF4 import Dataset
import glob
import time as tm
import csv
from collections import Counter
import scipy.ndimage as ndimage
import re
import h5py
import datetime as dt
import PF_functions as fa

#==========================================================
# Namelist
#==========================================================

# Directory and filename for input IMERG data
# Should obtain from create_FiT_input_files.py
datadirIM  = "/uufs/chpc.utah.edu/common/home/varble-group2/IMERG/"
fileidIM   = "3B-HHR.MS.MRG.3IMERG."
datahdf5   = True
datanc4    = False

# Directory and filename for FiT input files
# Should obtain from create_FiT_input_files.py
datadirFin = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_102519/FiT_input_files_2018_v2/"
fileidFin  = "IMERG_FiT_tholds_CapVer_"

# Directory and filename for FiT output data files
datadirFi  = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_102519/FiT_output_files_2018_v2/"
fileidFi1  = "IMERG_tracked_"
fileidFi2  = "_4Dobjects.nc"
fileidtxt  = "4Dobject_tree.txt"

# Directory and filename for output PF nc files
datadirout = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_102519/IMERGPFs_2018_v2/"
fileidout1 = "IPF_"
fileidout2 = ".nc"

# Only process a set of objects
subsetobs = True
obid1     = 100000
obid2     = 100100

# Specify merging distance from tracking (number of pixels)
mergedist = 30

# Specify reference time in format YYYY-MM-DD hh:mm:ss
reftime = "1900-01-01 00:00:00"

#==========================================================
# Initialize timer
#==========================================================

startnow = tm.time()

#==========================================================
# Get universal information
#==========================================================

# Open file and read information
print("Obtaining domain information")
datasetFin = Dataset(datadirFin+fileidFin+"00000.nc")
region = datasetFin.region
tholds = datasetFin.tholds
datestart = datasetFin.datestart
dateend = datasetFin.dateend
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

# Make list of filenames to search for and then find all
print("Generating list of files")

# Generate a list of filenames with dates to search for
start = dt.datetime.strptime(datestart, "%Y%m%d")
end = dt.datetime.strptime(dateend, "%Y%m%d")
datearr = (start + dt.timedelta(days=x) for x in range(0, (end-start).days))
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
lenIdir = len(datadirIM)

#==========================================================
# Read FiT text files and begin loop over object
#==========================================================

# Read text file
f       = open(datadirFi+fileidFi1+fileidtxt,'r') 
FiTinfo = np.genfromtxt(f, dtype="str",delimiter="\t",
                           skip_header=1)

# Create dictionary with number of occurrences of each IPF
obids = Counter([int(i) for i in FiTinfo[:,0]])
objs  = list(obids.keys())

# Subset by range of object ids
if subsetobs:
  objs = objs[obid1-100000:obid2-100000]

# Loop through all objects
for o in objs:
  print("Working on object "+str(o)+" with "+str(
           obids.get(o))+" times")

#==========================================================
# Begin loop over times for object
#==========================================================
  
  # Find all time indices and times relevant for object
  indob = [i for i, x in enumerate(list(
             map(int,FiTinfo[:,0]))) if x == int(o)]
  times = list(map(int,FiTinfo[indob,2]))

  # Loop over all times relevant to object
  for t in range(0,len(times)):
    print("Working on time: "+str(times[t]))

#==========================================================
# Read FiT objects file and process for specific object
#==========================================================

    # Obtain time from FiT input file
    Finfilename = datadirFin+fileidFin+str(
                    times[t]).zfill(5)+".nc"
    datasetFin  = Dataset(Finfilename)
    timenow     = str(datasetFin.time)

    # Open file and read map of object ids
    FiTfilename = datadirFi+fileidFi1+str(
                    times[t]).zfill(5)+fileidFi2
    print("Reading FiT object file: "+FiTfilename)
    datasetFi   = Dataset(FiTfilename)
    objects     = np.squeeze(
                    datasetFi.variables["value"][:,:,:])

    # Make obbin only 1 for specific object id
    obbin       = np.where(objects==int(o),1,0)

    # Find all location indices of object
    locindob    = np.where(obbin==1)

#==========================================================
# Read IMERG data and process for specific object
#==========================================================

    print("Reading IMERG data file: "+filenames[times[t]])
    
    if datanc4:
      # Open file
      datasetIM = Dataset(filenames[times[t]])

      # Read variables
      rain  = datasetIM.variables["precipitationCal"][:,:,:]

    if datahdf5:
      # Open file
      datasetIM = h5py.File(filenames[times[t]],'r')

      # Read variables
      rain = datasetIM["/Grid/precipitationCal"][:,:,:]

    # Change dimension order of variable
    raint = np.squeeze(np.transpose(rain,(0, 2, 1)))
    del(rain)
    rain = raint
    del(raint)

    # Subset data
    if ssreg:

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

      # Subset main dataset by indices
      latsub  = lat[idyS:idyN]
      lonsub  = lon[idxW:idxE]
      rainsub = rain[idyS:idyN,idxW:idxE]

      # Replace old vars with new subsetting vars
      del(rain,lat,lon)
      rain = rainsub
      lat = latsub
      lon = lonsub
      del(rainsub,latsub,lonsub)

    if sstd is not None:
      rainS = ndimage.gaussian_filter(rain,sigma=int(sstd))

    if swid is not None:
      rainS = ndimage.uniform_filter(rain,size=int(swid))

#==========================================================
# Calculate parameters and assign data to location
#==========================================================
    
    # Preallocate arrays for object f
    if t==0:
      datetime1       = [0] * len(times) # Times (YYYYMMDDhhmm)
      pieces1         = [0] * len(times) # No of pieces the object is made of
      centrallocx1    = [0] * len(times) # Central x-location of object in grid
      centrallocy1    = [0] * len(times) # Central y-location of object in grid
      centrallat1     = [0] * len(times) # Central latitude of object
      centrallon1     = [0] * len(times) # Central longitude of object
      time1           = [0.]* len(times) # Times in hours since reftime
 
      # Assign some simple values
      datetime1[t]    = int(timenow)
      pieces1[t]      = int(FiTinfo[indob,7][t].split(':')[1])
      centrallocx1[t] = int(re.split(':|,',FiTinfo[indob,6][t])[1])
      centrallocy1[t] = int(re.split(':|,',FiTinfo[indob,6][t])[2])
      centrallat1[t]  = lat[centrallocy1[t]]
      centrallon1[t]  = lon[centrallocx1[t]]

      # Make reference time object and calculate time since reftime
      currentdate = dt.datetime.strptime(timenow, "%Y%m%d%H%M")
      sincedate = dt.datetime.strptime(reftime, "%Y-%m-%d %H:%M:%S")
      delta = (currentdate - sincedate)
      time1[t] = delta.total_seconds()/3600.

      # Make dictionaries for data corresponding to PF at first time
      lats1   = {str(datetime1[t]): lat[locindob[0]]}
      lons1   = {str(datetime1[t]): lon[locindob[1]]}
      instrain1  = {str(datetime1[t]): rain[locindob[0],locindob[1]]}
      instrainS1 = {str(datetime1[t]): rainS[locindob[0],locindob[1]]}

    else:

      # Assign some simple values
      datetime1[t]    = int(timenow)
      pieces1[t]      = int(FiTinfo[indob,7][t].split(':')[1])
      centrallocx1[t] = int(re.split(':|,',FiTinfo[indob,6][t])[1])
      centrallocy1[t] = int(re.split(':|,',FiTinfo[indob,6][t])[2])
      centrallat1[t]  = lat[centrallocy1[t]]
      centrallon1[t]  = lon[centrallocx1[t]]

      # Calculate time in hours since reftime
      currentdate = dt.datetime.strptime(timenow, "%Y%m%d%H%M")
      delta = (currentdate - sincedate)
      time1[t] = delta.total_seconds()/3600.

      # For each new time, add to dictionary as a new key
      lats1[str(datetime1[t])]   = lat[locindob[0]]
      lons1[str(datetime1[t])]   = lon[locindob[1]]
      instrain1[str(datetime1[t])]  = rain[locindob[0],locindob[1]]
      instrainS1[str(datetime1[t])] = rainS[locindob[0],locindob[1]]

  # For PF, if all rain rates at all times are actually zero
  #  move onto the next PF
  for i in instrain1.keys():
    if len(instrain1[i][instrain1[i]>0])==0:
      a = True
    else:
      a = False
      break
  if a:
    print("All times in PF have zero precip, moving on")
    continue
  
  # For all other PFs, if all rain rates at any given time are
  #  zero, add that index to a list of indices to delete
  delind = []
  ind = 0
  for i in instrain1.keys():
    if len(instrain1[i][instrain1[i]>0])==0:
      delind.append(ind)
    ind = ind+1

  # Delete the key, value pairs in the dictionaries
  for i in delind:
    del(lats1[str(datetime1[i])])
    del(lons1[str(datetime1[i])]) 
    del(instrain1[str(datetime1[i])])
    del(instrainS1[str(datetime1[i])])

  # Delete the indices in the lists
  delind.reverse()
  for i in delind:
    del(datetime1[i])
    del(time1[i])
    del(pieces1[i])
    del(centrallocx1[i])
    del(centrallocy1[i])
    del(centrallat1[i])
    del(centrallon1[i])

#==========================================================
# Write NetCDF file for object
#==========================================================  

  # Creating netcdf file to write to
  filename = fileidout1+str(o)+"_"+str(datetime1[0])+"_"+str(
               int(round(centrallat1[0]))).zfill(2)+"_"+str(
               int(round(centrallon1[0]))).zfill(2)+".nc"

  print("Creating file: "+datadirout+filename)
  fileout = Dataset(datadirout+filename,'w',format='NETCDF4')

  # Global Attributes
  fileout.long_name   = 'Precipitation feature information'
  fileout.description = 'Information for a precipitation feature tracked by FiT algorithm with specifications as in the attributes below'
  fileout.history     = 'Created ' + tm.ctime(tm.time())
  fileout.source      = 'IMERG V06B: '+source
  fileout.objectID    = str(o)
  if ssreg:
    fileout.track_dom_latN   = latN
    fileout.track_dom_latS   = latS
    fileout.track_dom_lonW   = lonW
    fileout.track_dom_lonE   = lonE
    fileout.track_dom_region = region 
  if sstd is not None:
    fileout.smthcommand      = "ndimage.gaussian_filter(rain, sigma=int(smthstdd))"
    fileout.smthstdd         = stdd
  if swid is not None:
    fileout.smthcommand      = "ndimage.uniform_filter(rain, size=int(smthwidth))"
    fileout.smthwidth        = swid
    fileout.smthwidth_units  = "pixels"
    fileout.smth_description = "Smoothing used for tracking objects in FiT algorithm"
  fileout.splitdist   = mergedist
  fileout.splitdist_units = "pixels"
  fileout.splitdist_description = "Maximum distance objects can be separated but remain the same object after splitting in the FiT algorithm"
  fileout.tholds      = tholds
  fileout.tholds_units = "mm/hr"
  fileout.tholds_description = "Rain rate thresholds used for tracking with FiT algorithm"

  # Create dimensions in file
  time = fileout.createDimension('time', len(time1))

  # Create variables in file
  fa.write_var("time","Time","","time",np.float64,"hours since "+reftime,
              fileout,time1,datadirout+filename)

  fa.write_var("datetime","Date and time","","time",int,"YYYYMMDDhhmm",
              fileout,datetime1,datadirout+filename)

  description = "x-location of PF centroid in tracked domain"
  fa.write_var("centrallocx","Central x-location",description,"time",int,"",
              fileout,centrallocx1,datadirout+filename)

  description = "y-location of PF centroid in tracked domain"
  fa.write_var("centrallocy","Central y-location",description,"time",int,"",
              fileout,centrallocx1,datadirout+filename)

  fa.write_var("centrallat","Central latitude","Latitude of PF centroid","time",
                np.float64,"degreesNorth",fileout,centrallat1,datadirout+filename)

  fa.write_var("centrallon","Central longitude","Longitude of PF centroid","time",
                np.float64,"degreesEast",fileout,centrallon1,datadirout+filename)

  fa.write_var("pieces","Pieces","Number of pieces that make up the PF","time",
                int,"",fileout,pieces1,datadirout+filename)

  format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time."
  fa.write_group("lats","Latitudes","Latitudes of IMERG grid cell centers for PF",
                 "DegreesNorth",format1,fileout,lats1,datadirout+filename)

  fa.write_group("lons","Longitudes","Longitudes of IMERG grid cell centers for PF",
                 "DegreesEast",format1,fileout,lons1,datadirout+filename)

  description = "Instantaneous rain rates at IMERG grid cells corresponding to all latitude and longitude values in PF."
  fa.write_group("instrain","Instantaneous rain rate",description,
                 "mm/hr",format1,fileout,instrain1,datadirout+filename)

  if (sstd is not None) or (swid is not None):
    description = "Smoothed instantaneous rain rates at IMERG grid cells corresponding to all latitude and longitude values in PF"
    fa.write_group("instrainS","Smoothed instantaneous rain rate",description,
                 "mm/hr",format1,fileout,instrainS1,datadirout+filename)

  # Close file
  fileout.close()

#==========================================================
# End loop over objects
#==========================================================   
  
  del(indob,times,lats1,lons1,instrain1,instrainS1,datetime1,time1)

#==========================================================
# End timer
#==========================================================

endnow = tm.time()
print("Program took "+str(endnow - startnow)+"s")

#==========================================================
# End code
#==========================================================

