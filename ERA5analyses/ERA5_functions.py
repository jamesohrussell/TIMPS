#==================================================================
# Functions used to add variables to IMERG PFs
# Author: James Russell 2020
#
# * get_E5_subset_2D
#   - Gets a 2D subset of ERA5 data given central times, 
#      lats, and lons, for a specified variable.
#
# * get_E5_subset_file
#   - Given a time, finds the ERA5 file that is relevant, and
#      outputs the indices corresponding to the closest time 
#      coordinate within that file
#
# * get_E5_subset_2D_coords
#   - Gets the lat and lon indices and the coordinates within
#      a certain area surrounding a central lon and lat
#
# * get_E5_subset_2D_var
#   - Gets a subset of a specific ERA5 variable corresponding
#      to the time, latitude, and longitude coordinates provided
#
#==================================================================

#==================================================================
# Get subset of a variable from ERA5 files
#==================================================================

def get_E5_subset_2D(datadir,fileid,varname,timestr,clon,clat,hda):
  """
  Get a 3D subset of ERA5 data.

  Input: 
   1,2) Data directory anhd filename identifier
   3) A string of the variable name in ERA5 e.g."CAPE"
   4) A string indicating the current time in format:
       yyyy-mm-dd HH:MM:SS
   5,6) The central longitude and latitude for the subset
   7) Half the size of the subset in degrees (i.e. for a 10x10 
       degree subset hda==5)

  Output:
   1,2,3) Returns flattened lists of the ERA5 variable, and it's
           corresponding longitude and latitude coordinates
   

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/), scipy 1.2.1 (conda install 
   -c anaconda scipy; https://pypi.org/project/scipy/)
  """

  # Get the file handle and time coordinates
  fh,timi,times,ctime = get_E5_ss_file(
   datadir,fileid,timestr)

  # Find lat and lon coordinates and incices
  loni,lati,lonE5,latE5 = get_E5_ss_3D_coords(
   fh,clon,clat,hda)

  # Read in subset of specific variable
  varE5 = get_E5_ss_3D_var(
   fh,varname,timi,loni,lati,times,ctime)

  # Return data
  return(varE5,lonE5,latE5)



#==================================================================
# Get subset of a variable from ERA5 files
#==================================================================

def get_E5_ss_2D_fiti(datadir,fileid,timestr):
  """
  Find the file and output the time indices given a time

  Input: 
   1,2) Data directory anhd filename identifier
   3) A string indicating the current time in format:
       yyyy-mm-dd HH:MM:SS

  Output:
   1) A handle for the datafile
   2) Indices for the time coordinates
   3) A list of the times from the file
   4) The time to inerpolate to converted to ERA5 time units

  Requires 
  """

  # Import libraries
  import glob
  from netCDF4 import Dataset
  import TIPS_functions as fns

  # Select which file the time is within
  allfiles = sorted(glob.glob(datadir+fileid+"*"))
  ds0 = Dataset(allfiles[0])
  ctime = fns.time_since(timestr,ds0.variables["time"].units)
  for fi in allfiles:
    fh = Dataset(fi)
    if fh.variables["time"][0]<=ctime<=fh.variables["time"][-1]:
      break

  # Select the index(es) of the relevant time(s) 
  times = list(fh.variables["time"][:])
  if ctime in times:
    timi = times.index(ctime)
  else:
    timi = fns.k_closest(times,ctime,2)

  return(fh,timi,times,ctime)

#==================================================================
# Get subset of a variable from ERA5 files
#==================================================================

def get_E5_ss_2D_coords(fh,clon,clat,hda):
  """
  Find the horizontal coordinates and indices for a subset of 
   ERA5 data

  Input: 
   1) A handle for the datafile
   2,3) The central longitude and latitude for the subset
   4) Half the size of the subset in degrees (i.e. for a 10x10 
       degree subset hda==5)

  Output:
   1,2) Indices for the longitude and latitude coordinates
   3,4) Flattened lists of the longitude and latitude coordinates 

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/)
  """

  # Import libraries
  import numpy as np
  import TIPS_functions as fns

  # Find lat and lon indices
  if clon<0: clon = clon+360
  lat = fh.variables["latitude"][:]
  lon = fh.variables["longitude"][:]
  lati  = np.squeeze([fns.k_closest(lat,clat+hda,1), 
                      fns.k_closest(lat,clat-hda,1)])
  loni  = np.squeeze([fns.k_closest(lon,clon-hda,1),
                      fns.k_closest(lon,clon+hda,1)])
     
  # Get coordinates of that subset
  coords = np.meshgrid(
   fh.variables["latitude"][lati[0]:lati[1]+1],
   fh.variables["longitude"][loni[0]:loni[1]+1])

  # Return data
  return(loni,lati,coords[1].flatten(),coords[0].flatten())



#==================================================================
# Get subset of a variable from ERA5 files
#==================================================================

def get_E5_ss_2D_var(fh,varname,timi,loni,lati,times,ctime):
  """
  Get a 2D subset of ERA5 data.

  Input: 
   1) A file handle for the file
   2) A string of the variable name in ERA5 e.g."CAPE"
   3) The indices corresponding to time
   4,5) The indices corresponding to longitude and latitude ranges
   6,7) Times and the central time for interpolation

  Output:
   1) Returns a flattened list of the ERA5 variable corresponding 
    to coordinates from get_E5_subset_2D_coords

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/), scipy 1.2.1 (conda install 
   -c anaconda scipy; https://pypi.org/project/scipy/)
  """

  # Import libraries
  import numpy as np
  from scipy.interpolate import interp1d

  # Read in subset of data (interpolate in time if necessary)
  if hasattr(timi,"__len__"):
    varall = np.array(fh.variables[varname][timi[0]:timi[1]+1,
              lati[0]:lati[1]+1,loni[0]:loni[1]+1])
    varss  = interp1d([times[timi[0]],times[timi[1]]],
              varall,axis=0)(ctime)
  else:
    varss  = np.array(fh.variables[varname][timi,
              lati[0]:lati[1]+1,loni[0]:loni[1]+1])

  # Return data
  return(varss.flatten())

#==================================================================
# Get subset of a variable from ERA5 files
#==================================================================

def get_E5_ss_4D(datadir,fileid,varname,timestr1,timestr2,
                 lev1,lev2,lon1,lon2,lat1,lat2):
  """
  Get a 4D subset of ERA5 data.

  Input: 
   1)

  Output:
   1) 

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/), scipy 1.2.1 (conda install 
   -c anaconda scipy; https://pypi.org/project/scipy/)
  """

  files,timi1,timi2 = get_E5_ss_4D_fiti(datadir,fileid,
                                        timestr1,timestr2)

  

#==================================================================
# Get subset of a variable from ERA5 files
#==================================================================

def get_E5_ss_4D_fiti(datadir,fileid,timestr1,timestr2):
  """
  Find all files with data between two times and output the 
   indices of the times within the first and last files

  Input: 
   1,2) Data directory and file identifier
   3,4) Strings indicating the first and last times in format:
       yyyy-mm-dd HH:MM:SS

  Output:
   1) A list of strings indicating the paths and files
   2,3) Either a scalar or list of indices indicating the time 
    index(es) within that file corresponding to either the first
    or last times

  Requires glob and netCDF4
  """

  # Import libraries
  import glob
  from netCDF4 import Dataset
  import TIPS_functions as fns

  # Find all files 
  allfiles = sorted(glob.glob(datadir+fileid+"*"))

  # Convert first and last time to same units
  ds0 = Dataset(allfiles[0])
  time1 = fns.time_since(timestr1,ds0.variables["time"].units)
  time2 = fns.time_since(timestr2,ds0.variables["time"].units)

  # Loop over all files and initialize variables
  append = False
  ssfiles = []
  for fi in allfiles:
    fh = Dataset(fi)

    if append:
      times.extend(fh.variables["time"][:])

    # If first time is within current file set to append filenames
    if fh.variables["time"][0]<=time1<=fh.variables["time"][-1]:
      append = True

      # Select the index(es) of the relevant time(s)
      times1 = list(fh.variables["time"][:])
      if time1 in times1:
        timi1 = times1.index(time1)
      else:
        timi1 = fns.k_closest(times1,time1,1)[0]

      # Get times
      times = times1[timi1:]

    # Append files only if between first and last time
    if append:
      ssfiles.append(fi)

    # If last time is within current file
    if fh.variables["time"][0]<=time2<=fh.variables["time"][-1]:
    
      # Select the index(es) of the relevant time(s) 
      times2 = list(fh.variables["time"][:])
      if time2 in times2:
        timi2 = times2.index(time2)
      else:
        timi2 = fns.k_closest(times2,time2,1)[0]

      # Get times
      times.extend(times2[0:timi2+1])

      # Break the loop
      break

  # Return data
  return(ssfiles,timi1,timi2,times)



#==================================================================
# Get subset of a variable from ERA5 files
#==================================================================

def get_E5_ss_4D_coords(fh,lon1,lon2,lat1,lat2):
  """
  Find the horizontal coordinates and indices for a subset of 
   ERA5 data

  Input: 
   1) A handle for the datafile
   2,3) The central longitude and latitude for the subset
   4) Half the size of the subset in degrees (i.e. for a 10x10 
       degree subset hda==5)

  Output:
   1,2) Indices for the longitude and latitude coordinates
   3,4) Flattened lists of the longitude and latitude coordinates 

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/)
  """

  # Import libraries
  import numpy as np
  import TIPS_functions as fns

  # Find lat and lon indices
  if lon1<0: lon1 = lon1+360
  if lon2<0: lon2 = lon2+360
  print(lon1)
  print(lon2)
  lat = fh.variables["latitude"][:]
  lon = fh.variables["longitude"][:]
  lati  = np.squeeze([fns.k_closest(lat,lat2,1), 
                      fns.k_closest(lat,lat1,1)])
  loni  = np.squeeze([fns.k_closest(lon,lon1,1),
                      fns.k_closest(lon,lon2,1)])
  print(loni)

  # Get coordinates of that subset
  if loni[0]<loni[1]:
    coords = np.meshgrid(
     fh.variables["latitude"][lati[0]:lati[1]+1],
     fh.variables["longitude"][loni[0]:loni[1]+1])

  if loni[0]>loni[1]:
    coords = np.meshgrid(
     fh.variables["latitude"][lati[0]:lati[1]+1],
     list(fh.variables["longitude"][loni[0]:-1])+\
     list(fh.variables["longitude"][0:loni[1]+1]))

  # Return data
  return(loni,lati,coords[1].flatten(),coords[0].flatten())

