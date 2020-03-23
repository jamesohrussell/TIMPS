def driver_processFiTobs(o):
  "o corresponds to the object in objs"
  
  # Import python libraries (do not change)
  import numpy as np
  import h5py
  import csv
  import scipy.ndimage as ndimage
  import re
  from netCDF4 import Dataset
  import datetime as dt
  import time as tm
  import sys

#============================================================
# Read files output by main script
#============================================================

  # Read in namelist
  nl = Dataset("namelist_po.nc","r")

  # Load custom libraries
  sys.path.insert(0,nl.fnsdir)
  import misc_functions as mfns

  # Read text file
  f       = open(nl.datadirFi+nl.fileidFi1+nl.fileidtxt,'r') 
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

  # Loop over all times relevant to object
  for t in range(0,len(times)):

#============================================================
# Read FiT objects file and process for specific object
#============================================================

    # Open FiT input file
    Finfilename = nl.datadirFin+nl.fileidFin+str(
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
    FiTfilename = nl.datadirFi+nl.fileidFi1+str(
                    times[t]).zfill(5)+nl.fileidFi2
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
   
    if nl.datanc4=="True":
      # Open file
      datasetIM = Dataset(filenames[times[t]])

      # Read variables
      rain  = datasetIM.variables["precipitationCal"][:,:,:]

    if nl.datahdf5=="True":
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
      if nl.datanc4=="True":
        lat = datasetIM.variables["lat"][:]
        lon = datasetIM.variables["lon"][:]
      
      if nl.datahdf5=="True":
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

#============================================================
# Calculate parameters and assign data at first time
#============================================================
    
    # Preallocate arrays for object f
    if t==0:

      # Times (YYYYMMDDhhmm)
      datetime1       = [0] * len(times) 
      # No of pieces the object is made of
      pieces1         = [0] * len(times)
      # Central x-location of object in grid 
      centrallocx1    = [0] * len(times)
      # Central y-location of object in grid 
      centrallocy1    = [0] * len(times)
      # Central latitude of object
      centrallat1     = [0] * len(times) 
      # Central longitude of object
      centrallon1     = [0] * len(times)
      # Times in hours since reftime
      time1           = [0.]* len(times) 
 
      # Assign some simple values
      datetime1[t]    = int(timenow)
      pieces1[t]      = int(
       FiTinfo[indob,7][t].split(':')[1])
      centrallocx1[t] = int(re.split(':|,',
                            FiTinfo[indob,6][t])[1])
      centrallocy1[t] = int(re.split(':|,',
                            FiTinfo[indob,6][t])[2])
      centrallat1[t]  = lat[centrallocy1[t]]
      centrallon1[t]  = lon[centrallocx1[t]]

      # Make reference time object and calculate time
      #  in hours since reftime
      currentdate = dt.datetime.strptime(
       timenow,"%Y%m%d%H%M")
      sincedate = dt.datetime.strptime(
       nl.reftime, "%Y-%m-%d %H:%M:%S")
      delta = (currentdate - sincedate)
      time1[t] = delta.total_seconds()/3600.

      # Make dictionaries for data corresponding to PF
      #  at first time
      lats1   = {str(datetime1[t]): lat[locindob[0]]}
      lons1   = {str(datetime1[t]): lon[locindob[1]]}
      instrain1  = {str(datetime1[t]): rain[locindob[0],locindob[1]]}
      instrainS1 = {str(datetime1[t]): rainS[locindob[0],locindob[1]]}

#============================================================
# Calculate parameters and assign data at other times
#============================================================

    else:

      # Assign some simple values
      datetime1[t]    = int(timenow)
      pieces1[t]      = int(
       FiTinfo[indob,7][t].split(':')[1])
      centrallocx1[t] = int(re.split(
       ':|,',FiTinfo[indob,6][t])[1])
      centrallocy1[t] = int(re.split(
       ':|,',FiTinfo[indob,6][t])[2])
      centrallat1[t]  = lat[centrallocy1[t]]
      centrallon1[t]  = lon[centrallocx1[t]]

      # Calculate time in hours since reftime
      currentdate = dt.datetime.strptime(
       timenow,"%Y%m%d%H%M")
      delta = (currentdate - sincedate)
      time1[t] = delta.total_seconds()/3600.

      # For each new time, add to dictionary as a new key
      lats1[str(datetime1[t])] = lat[locindob[0]]
      lons1[str(datetime1[t])] = lon[locindob[1]]
      instrain1[str(datetime1[t])]  = \
       rain[locindob[0],locindob[1]]
      instrainS1[str(datetime1[t])] = \
       rainS[locindob[0],locindob[1]]

  # For PF, if all rain rates at all times are actually zero
  #  move onto the next PF
  for i in instrain1.keys():
    if len(instrain1[i][instrain1[i]>0])==0:
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

#============================================================
# Write NetCDF file for object
#============================================================

  # Creating netcdf file to write to
  filename = nl.fileidout1+str(objs[o])+"_"+ \
   str(datetime1[0])+"_"+ \
   str(int(round(centrallat1[0]))).zfill(2)+"_"+ \
   str(int(round(centrallon1[0]))).zfill(2)+".nc"

  print("Creating file: "+nl.datadirout+filename)
  fileout = Dataset(nl.datadirout+filename,
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
  fileout.splitdist   = int(nl.mergedist)
  fileout.splitdist_units = "pixels"
  fileout.splitdist_description = "Maximum distance objects can be separated but remain the same object after splitting in the FiT algorithm"
  fileout.tholds      = tholds
  fileout.tholds_units = "mm/hr"
  fileout.tholds_description = "Rain rate thresholds used for tracking with FiT algorithm"

  # Create dimensions in file
  time = fileout.createDimension('time', len(time1))

  # Create variables in file
  mfns.write_var("time","Time","","time",np.float64,
               "hours since "+nl.reftime,fileout,
               time1,nl.datadirout+filename)

  mfns.write_var("datetime","Date and time","","time",
               int,"YYYYMMDDhhmm",fileout,datetime1,
               nl.datadirout+filename)

#  description = "x-location of PF centroid in tracked domain"
#  mfns.write_var("centrallocx","Central x-location",
#               description,"time",int,"",fileout,
#               centrallocx1,nl.datadirout+filename)

#  description = "y-location of PF centroid in tracked domain"
#  mfns.write_var("centrallocy","Central y-location",
#               description,"time",int,"",fileout,
#               centrallocx1,nl.datadirout+filename)

  mfns.write_var("centrallat","Central latitude",
               "Latitude of PF centroid","time",
               np.float64,"degreesNorth",fileout,
               centrallat1,nl.datadirout+filename)

  mfns.write_var("centrallon","Central longitude",
               "Longitude of PF centroid","time",
               np.float64,"degreesEast",fileout,
               centrallon1,nl.datadirout+filename)

#  mfns.write_var("pieces","Pieces",
#               "Number of pieces that make up the PF",
#               "time",int,"",fileout,pieces1,
#               nl.datadirout+filename)

#  mfns.write_var("instrain","Instantaneous rain rate",
#                "IMERG instantaneous rain rates",
#                ("time","y","x"),np.float64,"mm/hr",
#                fileout,instrain1,nl.datadirout+filename)

  format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time."

  mfns.write_group("lats","Latitudes",
                 "Latitudes of IMERG grid cell centers for PF",
                 "DegreesNorth",format1,fileout,
                 lats1,nl.datadirout+filename)

  mfns.write_group("lons","Longitudes",
                 "Longitudes of IMERG grid cell centers for PF",
                 "DegreesEast",format1,fileout,
                 lons1,nl.datadirout+filename)

  description = "Instantaneous rain rates at IMERG grid cells corresponding to all latitude and longitude values in PF."
  mfns.write_group("instrain","Instantaneous rain rate",
                 description,"mm/hr",format1,fileout,
                 instrain1,nl.datadirout+filename)

#  if (sstd is not None) or (swid is not None):
#    description = "Smoothed instantaneous rain rates at IMERG grid cells corresponding to all latitude and longitude values in PF"
#    mfns.write_group("instrainS",
#                   "Smoothed instantaneous rain rate",
#                   description,"mm/hr",format1,fileout,
#                   instrainS1,nl.datadirout+filename)

  # Close file
  fileout.close()

