def driver_addvars(fn):
  """
  Driver to calculate various variables for a single PF.
  fn corresponds to the line number in filenames_av.txt.
  """

#==================================================================
# Begin script
#==================================================================

  # Import libraries
  import TIPS_functions as fns
  import ERA5_functions as E5fns
  from netCDF4 import Dataset
  import numpy as np

  # Read in namelist variables
  nl = Dataset("namelist_av.nc","r")

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

  # Get all times 
  datakeys = datalat.__dict__.keys()

#==================================================================
# Preallocate arrays
#==================================================================

  # Dictionaries for lats and lons
  lonsE5 = {}
  latsE5 = {}
  lats = {}
  lons = {}
  instrain = {}

#==================================================================
# Begin loop over times
#==================================================================

  c = 0
  for k in datakeys:

    # Get data for current time
    lats[k] = datalat.getncattr(k)
    lons[k] = datalon.getncattr(k)
    instrain[k] = datarain.getncattr(k)

    # Get locations of all pixels with nonzero precipitation
    latsnzk = lats[k][instrain[k]>0]
    lonsnzk = lons[k][instrain[k]>0]
    instrainnzk = instrain[k][instrain[k]>0]

#==================================================================
# Get appropriate ERA5 files and times
#==================================================================

    # Get a filename
    if nl.addCAPEE5=="True":
      fileid1 = nl.fileCAPEE5id
    elif nl.addTCWVE5=="True":
      fileid1 = nl.fileTCWVE5id
    elif nl.addCCTOE5=="True":
      fileid1 = nl.fileCCTOE5id
    elif nl.addSHRFE5=="True":
      fileid1 = nl.fileUSRFE5id
    elif nl.addSHRBE5=="True":
      fileid1 = nl.fileUSRBE5id
    elif nl.addSPHFE5=="True":
      fileid1 = nl.fileSPHFE5id
    elif nl.addSPHBE5=="True":
      fileid1 = nl.fileSPHBE5id

    # Set a time string
    timestr = str(k)[0:4]+"-"+str(k)[4:6]+"-"+str(k)[6:8]+\
              " "+str(k)[8:10]+":"+str(k)[10:12]+":00"

    # Find file and time indices
    fh,timi,times = E5fns.get_E5_ss_2D_fiti(
                     nl.dataE5dir,fileid1,timestr)[0:3]

#==================================================================
# Get coordinates and indices
#==================================================================

    # Find coordinates and indices
    loni,lati,lonsE5[k],latsE5[k] = \
     E5fns.get_E5_ss_2D_coords(
      fh,dataclon[c],dataclat[c],nl.hda)
    fh.close()

    # Create x and y coordinates
    if addctarea=="True":
      xE5 = np.linspace(-nl.hda,nl.hda,len(lonsE5[k]))
      yE5 = np.linspace(-nl.hda,nl.hda,len(latsE5[k]))

#==================================================================
# Assign ERA5 CAPE data
#==================================================================

    if nl.addCAPEE5=="True":

      # Preallocate array
      if c==0:
        CAPEE5 = np.zeros((len(datakeys),len(yE5),len(xE5)))

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileCAPEE5id,timestr)
      
      # Get a subset of the CAPE data
      CAPEE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fh,"CAPE",timi,loni,lati,times,ctime)
      CAPEE5units = fh.variables["CAPE"].units

      fh.close()

#==================================================================
# Assign ERA5 TCWV data
#==================================================================

    if nl.addTCWVE5=="True":

      # Preallocate array
      if c==0:
        TCWVE5 = np.zeros((len(datakeys),len(yE5),len(xE5)))

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileTCWVE5id,timestr)
      
      # Get a subset of the CAPE data
      TCWVE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fh,"TCWV",timi,loni,lati,times,ctime)
      TCWVE5units = fh.variables["TCWV"].units
      fh.close()

#==================================================================
# Assign ERA5 SPHF data
#==================================================================

    if nl.addSPHFE5=="True":

      # Preallocate array
      if c==0:
        SPHFE5 = np.zeros((len(datakeys),len(yE5),len(xE5)))

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileSPHFE5id,timestr)
      
      # Get a subset of the SPHF data
      SPHFE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fh,"Q",timi,loni,lati,times,ctime)
      SPHFE5units = fh.variables["Q"].units
      fh.close()

#==================================================================
# Assign ERA5 SPHB data
#==================================================================

    if nl.addSPHBE5=="True":

      # Preallocate array
      if c==0:
        SPHBE5 = np.zeros((len(datakeys),len(yE5),len(xE5)))

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileSPHBE5id,timestr)
      
      # Get a subset of the SPHF data
      SPHBE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fh,"Q",timi,loni,lati,times,ctime)
      SPHBE5units = fh.variables["Q"].units
      fh.close()

#==================================================================
# Assign ERA5 SHRF data
#==================================================================

    if nl.addSHRFE5=="True":

      # Preallocate array
      if c==0:
        USRFE5 = np.zeros((len(datakeys),len(yE5),len(xE5)))
        VSRFE5 = np.zeros((len(datakeys),len(yE5),len(xE5)))
      
      # Find file and time indices
      fhU,timiU,timesU,ctimeU = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileUSRFE5id,timestr)
      fhV,timiV,timesV,ctimeV = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileVSRFE5id,timestr)

      # Get a subset of the SPHF data
      USRFE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fhU,"USHR",timiU,loni,lati,timesU,ctimeU)
      VSRFE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fhV,"VSHR",timiV,loni,lati,timesV,ctimeV)
      SHRFE5units = fhU.variables["USHR"].units
      fhU.close()
      fhV.close()

#==================================================================
# Assign ERA5 SHRB data
#==================================================================

    if nl.addSHRBE5=="True":

      # Preallocate array
      if c==0:
        USRBE5 = np.zeros((len(datakeys),len(yE5),len(xE5)))
        VSRBE5 = np.zeros((len(datakeys),len(yE5),len(xE5)))
      
      # Find file and time indices
      fhU,timiU,timesU,ctimeU = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileUSRBE5id,timestr)
      fhV,timiV,timesV,ctimeV = E5fns.get_E5_ss_2D_fiti(
                       nl.dataE5dir,nl.fileVSRBE5id,timestr)

      # Get a subset of the SPHB data
      USRBE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fhU,"USHR",timiU,loni,lati,timesU,ctimeU)
      VSRBE5[c,:,:] = E5fns.get_E5_ss_2D_var(
                 fhV,"VSHR",timiV,loni,lati,timesV,ctimeV)
      SHRBE5units = fhU.variables["USHR"].units
      fhU.close()
      fhV.close()

#==================================================================
# Assign ERA5 TCRW data
#==================================================================

    if addctarearainchk=="True" or addctmeanrainchk=="True" or \
       addinarearainchk=="True" or add addinmeanrainchk=="True":  

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
      CAPEE5mean = np.mean(CAPEE5[c,:,:],axis=(1,2))

    if addSPHFE5=="True":
      SPHFE5mean = np.mean(SPHFE5[c,:,:],axis=(1,2))

    if addSPHBE5=="True":
      SPHBE5mean = np.mean(SPHBE5[c,:,:],axis=(1,2))

    if addSHRFE5=="True":
      USRFE5mean = np.mean(USRFE5[c,:,:],axis=(1,2))
      VSRFE5mean = np.mean(VSRFE5[c,:,:],axis=(1,2))

    if addSHRBE5=="True":
      USRBE5mean = np.mean(USRBE5[c,:,:],axis=(1,2))
      VSRBE5mean = np.mean(VSRBE5[c,:,:],axis=(1,2))

    if addCCTOE5=="True":
      CCTOE5mean = np.mean(CCTOE5[c,:,:],axis=(1,2))

#==================================================================
# Open file to write data
#==================================================================

  fileout = Dataset(f,'a')

#==================================================================
# Write coordinate data for environment information to file
#==================================================================

  if nl.addCAPEE5=="True" or nl.addTCWVE5=="True" or \
     nl.addCCTOE5=="True" or nl.addSHRFE5=="True" or \
     nl.addSHRBE5=="True" or nl.addSPHFE5=="True" or \
     nl.addSPHFE5=="True":

    format1 = "Data is in attribute and value pairs of the subgroup data. Attributes correspond to the date and time in YYYYMMDDhhmm format. Values of those attributes are lists of the data at that time. Data here corresponds to the location set by the equivalent attribute and value pairs in the lats and lons group."

    description = "longitudes corresponding to ERA5 data"
    fns.write_group("lonsE5","ERA5 longitudes",description,
                    "degreesE",format1,fileout,lonsE5,f)

    description = "latitudes corresponding to ERA5 data"
    fns.write_group("latsE5","ERA5 latitudes",description,
                    "degreesN",format1,fileout,latsE5,f)

    try: xE51 = fileout.createDimension('xE5',len(xE5))
    except: print("xE5 already defined")
    try: yE51 = fileout.createDimension('yE5',len(yE5))
    except: print("yE5 already defined")

    description = "Corresponds to ERA5 data centered on precipitation system centroid. Negative is west. Positive is east."
    fns.write_var("xE5","ERA5 zonal distance from centroid",
     description,("xE5"),np.float64,"degrees",fileout,xE5,f)
    description = "Corresponds to ERA5 data centered on precipitation system centroid. Negative is south. Positive is north."
    fns.write_var("yE5","Meridional distance from centroid",
     description,("yE5"),np.float64,"degrees",fileout,yE5,f)

#==================================================================
# Write CAPE information to file
#==================================================================

  if nl.addCAPEE5=="True":

    description = "Surface-based CAPE from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("CAPE_E5",
     "ERA5 Convective Available Potential Energy",description,
     ("time","yE5","xE5"),np.float64,CAPEE5units,fileout,CAPEE5,f)

#==================================================================
# Write TCWV information to file
#==================================================================

  if nl.addTCWVE5=="True":

    description = "Total column water vapor from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("TCWV_E5","ERA5 Total Column Water Vapor",
     description,("time","yE5","xE5"),np.float64,TCWVE5units,
     fileout,TCWVE5,f)

#==================================================================
# Write SPHF information to file
#==================================================================

  if nl.addSPHFE5=="True":

    description = "Specific humidity between 850-200 hPa from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("SPHU_850-200_E5",
     "ERA5 850-200 hPa mean specific humidity",description,
     ("time","yE5","xE5"),np.float64,SPHFE5units,fileout,SPHFE5,f)

#==================================================================
# Write SPHB information to file
#==================================================================

  if nl.addSPHBE5=="True":

    description = "Specific humidity between 1000-850 hPa from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("SPHU_1000-850_E5",
     "ERA5 1000-850 hPa mean specific humidity",description,
     ("time","yE5","xE5"),np.float64,SPHBE5units,fileout,SPHBE5,f)

#==================================================================
# Write SHRF information to file
#==================================================================

  if nl.addSHRFE5=="True":

    description = "Zonal wind shear between 850-200 hPa from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("USHR_850-200_E5",
     "ERA5 850-200 hPa zonal wind shear",description,
     ("time","yE5","xE5"),np.float64,SHRFE5units,fileout,USRFE5,f)

    description = "Meridional wind shear between 850-200 hPa from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("VSHR_850-200_E5",
     "ERA5 850-200 hPa meridional wind shear",description,
     ("time","yE5","xE5"),np.float64,SHRFE5units,fileout,VSRFE5,f)

#==================================================================
# Write SHRB information to file
#==================================================================

  if nl.addSHRBE5=="True":

    description = "Zonal wind shear between 1000-850 hPa from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("USHR_1000-850_E5",
     "ERA5 1000-850 hPa zonal wind shear",description,
     ("time","yE5","xE5"),np.float64,SHRBE5units,fileout,USRBE5,f)

    description = "Meridional wind shear between 1000-850 hPa from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("VSHR_1000-850_E5",
     "ERA5 1000-850 hPa meridional wind shear",description,
     ("time","yE5","xE5"),np.float64,SHRBE5units,fileout,VSRBE5,f)

#==================================================================
# Write CCTO information to file
#==================================================================

  if nl.addCCTOE5=="True":

    description = "Total cloud cover from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    fns.write_var("CCTO_E5","ERA5 Total Cloud Cover",
     description,("time","yE5","xE5"),np.float64,CCTOE5units,
     fileout,CCTOE5,f)

#==================================================================
# Close current file
#==================================================================

  fileout.close()

#==================================================================
# End processing of current PF
#==================================================================
