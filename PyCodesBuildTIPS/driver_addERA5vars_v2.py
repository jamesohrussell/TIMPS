def driver_addvars(fn):
  """
  Driver to calculate various variables for a single PF.
  fn corresponds to the line number in filenames_av.txt.
  """

#==================================================================
# Begin script
#==================================================================

  # Import libraries
  from netCDF4 import Dataset
  import numpy as np
  import datetime as dt
  import sys

  # Read in namelist variables
  nl = Dataset("namelist_av.nc","r")

  # Import custom libraries
  sys.path.insert(0,nl.fnsdir)
  import time_functions as tfns
  import misc_functions as mfns
  import ERA5_functions_v2 as E5fns

  # Read in filename for PF
  for i, row in enumerate(open("filenames_av.txt")):
    if i==fn:
      f = row[:-1]

  print("Working with file "+str(f))

  # Open file and assign data
  fd = Dataset(f)
  datalat  = fd.groups["lats"].groups["data"]
  datalon  = fd.groups["lons"].groups["data"]
  datarain = fd.groups["instrain"].groups["data"]
  dataclat = fd.variables["centrallat"][:]
  dataclon = fd.variables["centrallon"][:]
  datadtim = fd.variables["datetime"][:]

#==================================================================
# Make list of times and get appropriate files for current PS
#==================================================================

  # Define list of times in file
  timestrs = [str(it)[0:4]+"-"+str(it)[4:6]+"-"+str(it)[6:8]+
              " "+str(it)[8:10]+":"+str(it)[10:12]+":00" 
              for it in datadtim]

  # Define datetime objects for first and last times
  fto = dt.datetime(int(timestrs[0][0:4]),int(timestrs[0][5:7]),
   int(timestrs[0][8:10]),hour=int(timestrs[0][11:13]),
   minute=int(timestrs[0][14:16]),second=int(timestrs[0][17:19]))
  lto = dt.datetime(int(timestrs[-1][0:4]),int(timestrs[-1][5:7]),
   int(timestrs[-1][8:10]),hour=int(timestrs[-1][11:13]),
   minute=int(timestrs[-1][14:16]),second=int(timestrs[-1][17:19]))
  
  # Find times before and after and add to time list
  timestrs = [str(fto-dt.timedelta(hours=int(it))) 
              for it in nl.hoursbefore] + timestrs + \
             [str(lto+dt.timedelta(hours=int(it))) 
                for it in nl.hoursafter]

  # Find range of files and times in those files
  files = E5fns.get_E5_ss_4D_fiti(
   nl.dataE5dir,nl.fileTCWVE5id,timestrs[0],timestrs[-1])[0]
  times = [tfns.time_since(i,fd.variables["time"].units)
            for i in timestrs]

#==================================================================
# Get corresponding lists of central latitude and longitude
#==================================================================

  dataclat = [round(i,2) for i in 
             [dataclat[0]]*len(nl.hoursbefore) + \
             [i for i in dataclat] + \
             [dataclat[-1]]*len(nl.hoursafter)]
  dataclon = [round(i,2) for i in 
             [dataclon[0]]*len(nl.hoursbefore) + \
             [i for i in dataclon] + \
             [dataclon[-1]]*len(nl.hoursafter)]

#==================================================================
# Begin loop over times
#==================================================================

  c = 0
  for k in timestrs:

#==================================================================
# Get appropriate ERA5 files and times
#==================================================================

    # Find file and time indices
    fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                     nl.dataE5dir,nl.fileTCWVE5id,k)

#==================================================================
# Get coordinates and indices
#==================================================================

    # Find coordinates and indices
    loni,lati,lons,lats = E5fns.get_E5_ss_2D_coords(
     fh,dataclat,dataclon,nl.hda)[0:2]
 
    # Close file
    fh.close()

    # Define x coordinate
    xE5 = np.linspace(-nl.hda,nl.hda,len(lonsE5[k]))
    yE5 = np.linspace(-nl.hda,nl.hda,len(latsE5[k]))

#==================================================================
# Assign ERA5 TCWV data
#==================================================================

    if nl.addTCWVE5=="True":

      # Preallocate array
      if c==0:
        TCWVE5 = np.zeros((len(times),len(yE5),len(xE5)))

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileTCWVE5id,timestr)
      
      # Get a subset of the CAPE data
      TCWVE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fh,"TCWV",timi,loni,lati,times,ctime)
      TCWVE5units = fh.variables["TCWV"].units
      fh.close()

#==================================================================
# Assign ERA5 TCRW data
#==================================================================

    if addctmeanrainchk=="True" or addinmeanrainchk=="True":  

      # Preallocate array
      if c==0:
        TCRWE5 = np.zeros((len(datakeys),len(yE5),len(xE5)))

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileTCRWE5id,timestr)
      
      # Get a subset of the TCRW data
      TCRWE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fh,"TCRW",timi,loni,lati,times,ctime)
      TCRWE5units = fh.variables["TCRW"].units
      fh.close()
      
#==================================================================
# End loops over objects and times
#==================================================================
      
    # Advance counter
    c = c + 1

  fd.close()

#==================================================================
# Calculate mean time-series
#==================================================================
      
  if addctmean=="True":

    if addCAPEE5=="True":
      TCWVE5mean = np.mean(TCWVE5[c,:,:],axis=(1,2))

#==================================================================
# Open file to write data
#==================================================================

  fileout = Dataset(f,'a')

#==================================================================
# Write coordinate data for environment information to file
#==================================================================

  if nl.addTCWVE5=="True":

    format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."

    description = "longitudes corresponding to ERA5 data"
    mfns.write_group("lonsE5","ERA5 longitudes",description,
                    "degreesE",format1,fileout,lonsE5,f)

    description = "latitudes corresponding to ERA5 data"
    mfns.write_group("latsE5","ERA5 latitudes",description,
                    "degreesN",format1,fileout,latsE5,f)

    try: xE51 = fileout.createDimension('xE5',len(xE5))
    except: print("xE5 already defined")
    try: yE51 = fileout.createDimension('yE5',len(yE5))
    except: print("yE5 already defined")

    description = "Corresponds to ERA5 data centered on precipitation system centroid. Negative is west. Positive is east."
    mfns.write_var("xE5","ERA5 zonal distance from centroid",
     description,("xE5"),np.float64,"degrees",fileout,xE5,f)
    description = "Corresponds to ERA5 data centered on precipitation system centroid. Negative is south. Positive is north."
    mfns.write_var("yE5","Meridional distance from centroid",
     description,("yE5"),np.float64,"degrees",fileout,yE5,f)

#==================================================================
# Write TCWV information to file
#==================================================================

  if nl.addTCWVE5=="True":

    description = "Total column water vapor from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("TCWV_E5","ERA5 Total Column Water Vapor",
     description,("time","yE5","xE5"),np.float64,TCWVE5units,
     fileout,TCWVE5,f)

#==================================================================
# Close current file
#==================================================================

  fileout.close()

#==================================================================
# End processing of current PF
#==================================================================
