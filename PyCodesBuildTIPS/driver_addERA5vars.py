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
  nl = Dataset("namelist_av_E5.nc","r")

  # Import custom libraries
  sys.path.insert(0,nl.fnsdir)
  import time_functions as tfns
  import misc_functions as mfns
  import ERA5_functions as E5fns

  # Read in filename for PF
  for i, row in enumerate(open("filenames_av_E5.txt")):
    if i==fn:
      f = row[:-1]

  print("Working with file "+str(f))

  # Open file and assign data
  fd = Dataset(f)
  dataclat = fd.variables["centrallat"][:]
  dataclon = fd.variables["centrallon"][:]
  datadtim = fd.variables["datetime"][:]
  timeunits = fd.variables["time"].units

  # Preallocate arrays
  lonsE5 = {}; latsE5 = {}
  files = {}

#==================================================================
# Make list of times and get list files for current TIPS
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
  tE5 = [tfns.time_since(i,timeunits) for i in timestrs]

  # Find files
  if nl.addrainchk=="True": files["TCRW"] = E5fns.get_E5_ss_files(
    nl.dataE5dir,nl.fileTCRWE5id,timestrs[0],timestrs[-1])
  if nl.addTCWVE5=="True": files["TCWV"] = E5fns.get_E5_ss_files(
    nl.dataE5dir,nl.fileTCWVE5id,timestrs[0],timestrs[-1])
  if nl.addCAPEE5=="True": files["CAPE"] = E5fns.get_E5_ss_files(
    nl.dataE5dir,nl.fileCAPEE5id,timestrs[0],timestrs[-1])
  if nl.addSR18E5=="True": 
    files["USR18"] = E5fns.get_E5_ss_files(
     nl.dataE5dir,nl.fileUSR18E5id,timestrs[0],timestrs[-1])
    files["VSR18"] = E5fns.get_E5_ss_files(
     nl.dataE5dir,nl.fileVSR18E5id,timestrs[0],timestrs[-1])
  if nl.addSR82E5=="True": 
    files["USR82"] = E5fns.get_E5_ss_files(
     nl.dataE5dir,nl.fileUSR82E5id,timestrs[0],timestrs[-1])
    files["VSR82"] = E5fns.get_E5_ss_files(
     nl.dataE5dir,nl.fileVSR82E5id,timestrs[0],timestrs[-1])
  if nl.addSH18E5=="True": files["SH18"] = E5fns.get_E5_ss_files(
    nl.dataE5dir,nl.fileSH18E5id,timestrs[0],timestrs[-1])
  if nl.addSH82E5=="True": files["SH82"] = E5fns.get_E5_ss_files(
    nl.dataE5dir,nl.fileSH82E5id,timestrs[0],timestrs[-1])

#==================================================================
# Get coordinates and their indices
#==================================================================

  # Get corresponding lists of central latitude and longitude
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
# Get coordinates
#==================================================================

    # Case: Add centered area or centered mean with rain check
    if nl.addctarea=="True" or \
      (nl.addctmean=="True" and nl.addrainchk=="True"):
  
      # Find coordinates and indices
      keys = [str(k) for k in files.keys()]
      loni,lati,lonsE5[k],latsE5[k] = E5fns.get_E5_ss_2D_coords(
        Dataset(files[keys[0]][0]),dataclon[c],dataclat[c],nl.hda)

      # Assign new coords
      xE5 = np.linspace(-nl.hda,nl.hda,len(lonsE5[k]))
      yE5 = np.linspace(-nl.hda,nl.hda,len(latsE5[k]))

    # Case: Only add centered mean
    elif nl.addctmean=="True" and nl.addnorainchk=="True":
  
      # Find indices
      keys = [str(k) for k in files.keys()]
      loni,lati = E5fns.get_E5_ss_coords(
        Dataset(files[keys[0]][0]),dataclon[c],dataclat[c],nl.hda)

#==================================================================
# Preallocate arrays
#==================================================================

    if c==0:

      if nl.addTCWVE5=="True":

        if nl.addctarea=="True" and nl.addnorainchk=="True":
          TCWV_area_E5 = np.zeros((
           len(timestrs),len(yE5),len(xE5)))  

        if nl.addctmean=="True" and nl.addnorainchk=="True": 
          TCWV_mean_E5 = [-999]*len(timestrs)
  
        if nl.addctarea=="True" and nl.addrainchk=="True": 
          TCWV_area_nr_E5 = np.zeros((
           len(timestrs),len(yE5),len(xE5)))  

        if nl.addctmean=="True" and nl.addrainchk=="True": 
          TCWV_mean_nr_E5 = [-999]*len(timestrs)

      if nl.addCAPEE5=="True":

        if nl.addctarea=="True" and nl.addnorainchk=="True":
          CAPE_area_E5 = np.zeros((
           len(timestrs),len(yE5),len(xE5)))  

        if nl.addctmean=="True" and nl.addnorainchk=="True": 
          CAPE_mean_E5 = [-999]*len(timestrs)
  
        if nl.addctarea=="True" and nl.addrainchk=="True": 
          CAPE_area_nr_E5 = np.zeros((
           len(timestrs),len(yE5),len(xE5)))  

        if nl.addctmean=="True" and nl.addrainchk=="True": 
          CAPE_mean_nr_E5 = [-999]*len(timestrs)

      if nl.addSR18E5=="True":

        if nl.addctarea=="True" and nl.addnorainchk=="True":
          USR18_area_E5 = np.zeros((
           len(timestrs),len(yE5),len(xE5))) 
          VSR18_area_E5 = np.zeros((
           len(timestrs),len(yE5),len(xE5)))

        if nl.addctmean=="True" and nl.addnorainchk=="True": 
          USR18_mean_E5 = [-999]*len(timestrs)
          VSR18_mean_E5 = [-999]*len(timestrs)
  
        if nl.addctarea=="True" and nl.addrainchk=="True": 
          USR18_area_nr_E5 = np.zeros((
           len(timestrs),len(yE5),len(xE5)))  
          VSR18_area_nr_E5 = np.zeros((
           len(timestrs),len(yE5),len(xE5)))  

        if nl.addctmean=="True" and nl.addrainchk=="True": 
          USR18_mean_nr_E5 = [-999]*len(timestrs)
          VSR18_mean_nr_E5 = [-999]*len(timestrs)

      if nl.addSR82E5=="True":

        if nl.addctarea=="True" and nl.addnorainchk=="True":
          USR82_area_E5 = np.zeros((
           len(timestrs),len(yE5),len(xE5))) 
          VSR82_area_E5 = np.zeros((
           len(timestrs),len(yE5),len(xE5)))

        if nl.addctmean=="True" and nl.addnorainchk=="True": 
          USR82_mean_E5 = [-999]*len(timestrs)
          VSR82_mean_E5 = [-999]*len(timestrs)
  
        if nl.addctarea=="True" and nl.addrainchk=="True": 
          USR82_area_nr_E5 = np.zeros((
           len(timestrs),len(yE5),len(xE5)))  
          VSR82_area_nr_E5 = np.zeros((
           len(timestrs),len(yE5),len(xE5)))  

        if nl.addctmean=="True" and nl.addrainchk=="True": 
          USR82_mean_nr_E5 = [-999]*len(timestrs)
          VSR82_mean_nr_E5 = [-999]*len(timestrs)

      if nl.addSH18E5=="True":

        if nl.addctarea=="True" and nl.addnorainchk=="True":
          SH18_area_E5 = np.zeros((
           len(timestrs),len(yE5),len(xE5)))  

        if nl.addctmean=="True" and nl.addnorainchk=="True": 
          SH18_mean_E5 = [-999]*len(timestrs)
  
        if nl.addctarea=="True" and nl.addrainchk=="True": 
          SH18_area_nr_E5 = np.zeros((
           len(timestrs),len(yE5),len(xE5)))  

        if nl.addctmean=="True" and nl.addrainchk=="True": 
          SH18_mean_nr_E5 = [-999]*len(timestrs)

      if nl.addSH82E5=="True":

        if nl.addctarea=="True" and nl.addnorainchk=="True":
          SH82_area_E5 = np.zeros((
           len(timestrs),len(yE5),len(xE5)))  

        if nl.addctmean=="True" and nl.addnorainchk=="True": 
          SH82_mean_E5 = [-999]*len(timestrs)
  
        if nl.addctarea=="True" and nl.addrainchk=="True": 
          SH82_area_nr_E5 = np.zeros((
           len(timestrs),len(yE5),len(xE5)))  

        if nl.addctmean=="True" and nl.addrainchk=="True": 
          SH82_mean_nr_E5 = [-999]*len(timestrs)

      if nl.addrainchk=="True":
        TCRW_area_E5 = np.zeros((
         len(timestrs),len(yE5),len(xE5)))  

#==================================================================
# Get rain check data
#==================================================================

    if nl.addrainchk=="True":
      # Find file and time indices
      fhR,timiR,timesR,ctimeR = E5fns.get_E5_ss_2D_fiti(
       files["TCRW"],timestrs[c])

      TCRWE5 = E5fns.get_E5_ss_2D_var(
       fhR,"TCRW",timiR,loni,lati,timesR,ctimeR)

      if nl.addctarea=="True":
        TCRW_area_E5[c,:,:] = TCRWE5
        TCRWE5units = fhR[0].variables["TCRW"].units

      # Close file
      for fii in fhR:
        fii.close()

#==================================================================
# Get TCWV data
#==================================================================

    if nl.addTCWVE5=="True":

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       files["TCWV"],timestrs[c])

      # Get units
      TCWVE5units = fh[0].variables["TCWV"].units

      # Get an area
      if nl.addctarea=="True" or nl.addctmean=="True":
        TCWVE5 = E5fns.get_E5_ss_2D_var(
                 fh,"TCWV",timi,loni,lati,times,ctime)

      # Assign to area centered on PF
      if nl.addctarea=="True" and nl.addnorainchk=="True":
        TCWV_area_E5[c,:,:] = TCWVE5

      # Mean of area centered on PF
      if nl.addctmean=="True" and nl.addnorainchk=="True":
        TCWV_mean_E5[c] = np.mean(TCWVE5)

      # Calculate without raining pixels
      if (nl.addctarea=="True" or nl.addctmean=="True") \
        and nl.addrainchk=="True": 
        TCWVE5_nr = np.where(TCRWE5>0.001,np.nan,TCWVE5)

      # Area centered on PF without rain
      if nl.addctarea=="True" and nl.addrainchk=="True":
        TCWV_area_nr_E5[c,:,:] = TCWVE5_nr

      # Mean of area centered on PF
      if nl.addctmean=="True" and nl.addrainchk=="True":
        varnan = TCWVE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<nl.avgmissfrac:
          TCWV_mean_nr_E5[c] = np.nanmean(TCWVE5_nr)
      
      # Close file
      for fii in fh:
        fii.close()

#==================================================================
# Get CAPE data
#==================================================================

    if nl.addCAPEE5=="True":

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       files["CAPE"],timestrs[c])

      # Get units
      CAPEE5units = fh[0].variables["CAPE"].units

      # Get an area
      if nl.addctarea=="True" or nl.addctmean=="True":
        CAPEE5 = E5fns.get_E5_ss_2D_var(
                 fh,"CAPE",timi,loni,lati,times,ctime)

      # Assign to area centered on PF
      if nl.addctarea=="True" and nl.addnorainchk=="True":
        CAPE_area_E5[c,:,:] = CAPEE5

      # Mean of area centered on PF
      if nl.addctmean=="True" and nl.addnorainchk=="True":
        CAPE_mean_E5[c] = np.mean(CAPEE5)

      # Calculate without raining pixels
      if (nl.addctarea=="True" or nl.addctmean=="True") \
        and nl.addrainchk=="True": 
        CAPEE5_nr = np.where(TCRWE5>0.001,np.nan,CAPEE5)

      # Area centered on PF without rain
      if nl.addctarea=="True" and nl.addrainchk=="True":
        CAPE_area_nr_E5[c,:,:] = CAPEE5_nr

      # Mean of area centered on PF
      if nl.addctmean=="True" and nl.addrainchk=="True":
        varnan = CAPEE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<nl.avgmissfrac:
          CAPE_mean_nr_E5[c] = np.nanmean(CAPEE5_nr)
      
      # Close file
      for fii in fh:
        fii.close()

#==================================================================
# Get USR18 data
#==================================================================

    if nl.addSR18E5=="True":

      # Find file and time indices
      fhU,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       files["USR18"],timestrs[c])
      fhV,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       files["VSR18"],timestrs[c])

      # Get units
      SR18E5units = fhU[0].variables["USHR"].units

      # Get an area
      if nl.addctarea=="True" or nl.addctmean=="True":
        USR18E5 = E5fns.get_E5_ss_2D_var(
                 fhU,"USHR",timi,loni,lati,times,ctime)
        VSR18E5 = E5fns.get_E5_ss_2D_var(
                 fhV,"VSHR",timi,loni,lati,times,ctime)

      # Assign to area centered on PF
      if nl.addctarea=="True" and nl.addnorainchk=="True":
        USR18_area_E5[c,:,:] = USR18E5
        VSR18_area_E5[c,:,:] = VSR18E5

      # Mean of area centered on PF
      if nl.addctmean=="True" and nl.addnorainchk=="True":
        USR18_mean_E5[c] = np.mean(USR18E5)
        VSR18_mean_E5[c] = np.mean(VSR18E5)

      # Calculate without raining pixels
      if (nl.addctarea=="True" or nl.addctmean=="True") \
        and nl.addrainchk=="True":
        USR18E5_nr = np.where(TCRWE5>0.001,np.nan,USR18E5)
        VSR18E5_nr = np.where(TCRWE5>0.001,np.nan,VSR18E5)

      # Area centered on PF without rain
      if nl.addctarea=="True" and nl.addrainchk=="True":
        USR18_area_nr_E5[c,:,:] = USR18E5_nr
        VSR18_area_nr_E5[c,:,:] = VSR18E5_nr

      # Mean of area centered on PF
      if nl.addctmean=="True" and nl.addrainchk=="True":
        varnan = USR18E5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<nl.avgmissfrac:
          USR18_mean_nr_E5[c] = np.nanmean(USR18E5_nr)
        varnan = VSR18E5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<nl.avgmissfrac:
          VSR18_mean_nr_E5[c] = np.nanmean(VSR18E5_nr)
      
      # Close file
      for fii in fhU:
        fii.close()
      for fii in fhV:
        fii.close()

#==================================================================
# Get USR82 data
#==================================================================

    if nl.addSR82E5=="True":

      # Find file and time indices
      fhU,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       files["USR82"],timestrs[c])
      fhV,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       files["VSR82"],timestrs[c])

      # Get units
      SR82E5units = fhU[0].variables["USHR"].units

      # Get an area
      if nl.addctarea=="True" or nl.addctmean=="True":
        USR82E5 = E5fns.get_E5_ss_2D_var(
                 fhU,"USHR",timi,loni,lati,times,ctime)
        VSR82E5 = E5fns.get_E5_ss_2D_var(
                 fhV,"VSHR",timi,loni,lati,times,ctime)

      # Assign to area centered on PF
      if nl.addctarea=="True" and nl.addnorainchk=="True":
        USR82_area_E5[c,:,:] = USR82E5
        VSR82_area_E5[c,:,:] = VSR82E5

      # Mean of area centered on PF
      if nl.addctmean=="True" and nl.addnorainchk=="True":
        USR82_mean_E5[c] = np.mean(USR82E5)
        VSR82_mean_E5[c] = np.mean(VSR82E5)

      # Calculate without raining pixels
      if (nl.addctarea=="True" or nl.addctmean=="True") \
        and nl.addrainchk=="True":
        USR82E5_nr = np.where(TCRWE5>0.001,np.nan,USR82E5)
        VSR82E5_nr = np.where(TCRWE5>0.001,np.nan,VSR82E5)

      # Area centered on PF without rain
      if nl.addctarea=="True" and nl.addrainchk=="True":
        USR82_area_nr_E5[c,:,:] = USR82E5_nr
        VSR82_area_nr_E5[c,:,:] = VSR82E5_nr

      # Mean of area centered on PF
      if nl.addctmean=="True" and nl.addrainchk=="True":
        varnan = USR82E5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<nl.avgmissfrac:
          USR82_mean_nr_E5[c] = np.nanmean(USR82E5_nr)
        varnan = VSR82E5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<nl.avgmissfrac:
          VSR82_mean_nr_E5[c] = np.nanmean(VSR82E5_nr)
      
      # Close file
      for fii in fhU:
        fii.close()
      for fii in fhV:
        fii.close()

#==================================================================
# Get SH18 data
#==================================================================

    if nl.addSH18E5=="True":

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       files["SH18"],timestrs[c])

      # Get units
      SH18E5units = fh[0].variables["Q"].units

      # Get an area
      if nl.addctarea=="True" or nl.addctmean=="True":
        SH18E5 = E5fns.get_E5_ss_2D_var(
                 fh,"Q",timi,loni,lati,times,ctime)

      # Assign to area centered on PF
      if nl.addctarea=="True" and nl.addnorainchk=="True":
        SH18_area_E5[c,:,:] = SH18E5

      # Mean of area centered on PF
      if nl.addctmean=="True" and nl.addnorainchk=="True":
        SH18_mean_E5[c] = np.mean(SH18E5)

      # Calculate without raining pixels
      if (nl.addctarea=="True" or nl.addctmean=="True") \
        and nl.addrainchk=="True": 
        SH18E5_nr = np.where(TCRWE5>0.001,np.nan,SH18E5)

      # Area centered on PF without rain
      if nl.addctarea=="True" and nl.addrainchk=="True":
        SH18_area_nr_E5[c,:,:] = SH18E5_nr

      # Mean of area centered on PF
      if nl.addctmean=="True" and nl.addrainchk=="True":
        varnan = SH18E5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<nl.avgmissfrac:
          SH18_mean_nr_E5[c] = np.nanmean(SH18E5_nr)
      
      # Close file
      for fii in fh:
        fii.close()
     
#==================================================================
# Get SH82 data
#==================================================================

    if nl.addSH82E5=="True":

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       files["SH82"],timestrs[c])

      # Get units
      SH82E5units = fh[0].variables["Q"].units

      # Get an area
      if nl.addctarea=="True" or nl.addctmean=="True":
        SH82E5 = E5fns.get_E5_ss_2D_var(
                 fh,"Q",timi,loni,lati,times,ctime)

      # Assign to area centered on PF
      if nl.addctarea=="True" and nl.addnorainchk=="True":
        SH82_area_E5[c,:,:] = SH82E5

      # Mean of area centered on PF
      if nl.addctmean=="True" and nl.addnorainchk=="True":
        SH82_mean_E5[c] = np.mean(SH82E5)

      # Calculate without raining pixels
      if (nl.addctarea=="True" or nl.addctmean=="True") \
        and nl.addrainchk=="True": 
        SH82E5_nr = np.where(TCRWE5>0.001,np.nan,SH82E5)

      # Area centered on PF without rain
      if nl.addctarea=="True" and nl.addrainchk=="True":
        SH82_area_nr_E5[c,:,:] = SH82E5_nr

      # Mean of area centered on PF
      if nl.addctmean=="True" and nl.addrainchk=="True":
        varnan = SH82E5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<nl.avgmissfrac:
          SH82_mean_nr_E5[c] = np.nanmean(SH82E5_nr)
      
      # Close file
      for fii in fh:
        fii.close()

#==================================================================
# End loops over objects and times
#==================================================================
      
    # Advance counter
    c = c + 1

  fd.close()

#==================================================================
# Open file to write data
#==================================================================

  fileout = Dataset(f,'a')

#==================================================================
# Write coordinate data for environment information to file
#==================================================================

  try: tE51 = fileout.createDimension('tE5',len(tE5))
  except: print("tE5 already defined")

  description = "Corresponds to ERA5 data. Different from time dimension since times are added before and after the TIPS exists."
  mfns.write_var("tE5","ERA5 time",description,("tE5"),
    np.float64,timeunits,fileout,tE5,f,float(-999))

  if nl.addctarea=="True":

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
      description,("xE5"),np.float64,"degrees",fileout,xE5,f,
      float(-999))

    description = "Corresponds to ERA5 data centered on precipitation system centroid. Negative is south. Positive is north."
    mfns.write_var("yE5","Meridional distance from centroid",
      description,("yE5"),np.float64,"degrees",fileout,yE5,f,
      float(-999))

    if nl.addrainchk=="True" and nl.addctarea=="True":

      description = "Total column rain water from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("TCRW_area_E5","ERA5 Total Column Rain Water",
       description,("tE5","yE5","xE5"),np.float64,TCRWE5units,
       fileout,TCRW_area_E5,f,float(-999))

#==================================================================
# Write TCWV information to file
#==================================================================

  if nl.addTCWVE5=="True":

    if nl.addctarea=="True" and nl.addnorainchk=="True":
      description = "Total column water vapor from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("TCWV_area_E5","ERA5 Total Column Water Vapor",
       description,("tE5","yE5","xE5"),np.float64,TCWVE5units,
       fileout,TCWV_area_E5,f,float(-999))

    if nl.addctmean=="True"and nl.addnorainchk=="True":
      description = "Total column water vapor from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("TCWV_mean_E5",
       "ERA5 Mean Total Column Water Vapor",description,("tE5"),
       np.float64,TCWVE5units,fileout,TCWV_mean_E5,f,float(-999))

    if nl.addctarea=="True"and nl.addrainchk=="True":
      TCWV_area_nr_E5 = np.where(np.isnan(TCWV_area_nr_E5),float(-999),TCWV_area_nr_E5)
      description = "Total column water vapor from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("TCWV_area_nr_E5",
       "ERA5 Rain-Checked Total Column Water Vapor",
       description,("tE5","yE5","xE5"),np.float64,TCWVE5units,
       fileout,TCWV_area_nr_E5,f,float(-999))

    if nl.addctmean=="True"and nl.addrainchk=="True":
      description = "Total column water vapor from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("TCWV_mean_nr_E5",
       "ERA5 Rain-Checked Mean Total Column Water Vapor",
       description,("tE5"),np.float64,TCWVE5units,
       fileout,TCWV_mean_nr_E5,f,float(-999))

#==================================================================
# Write CAPE information to file
#==================================================================

  if nl.addCAPEE5=="True":

    if nl.addctarea=="True" and nl.addnorainchk=="True":
      description = "Convective Available Potential Energy from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("CAPE_area_E5","ERA5 CAPE",
       description,("tE5","yE5","xE5"),np.float64,CAPEE5units,
       fileout,CAPE_area_E5,f,float(-999))

    if nl.addctmean=="True"and nl.addnorainchk=="True":
      description = "Convective Available Potential Energy from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("CAPE_mean_E5","ERA5 Mean CAPE",
       description,("tE5"),np.float64,CAPEE5units,fileout,
       CAPE_mean_E5,f,float(-999))

    if nl.addctarea=="True"and nl.addrainchk=="True":
      CAPE_area_nr_E5 = np.where(np.isnan(CAPE_area_nr_E5),float(-999),CAPE_area_nr_E5)
      description = "Convective Available Potential Energy from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("CAPE_area_nr_E5","ERA5 Rain-Checked CAPE",
       description,("tE5","yE5","xE5"),np.float64,CAPEE5units,
       fileout,CAPE_area_nr_E5,f,float(-999))

    if nl.addctmean=="True"and nl.addrainchk=="True":
      description = "Convective Available Potential Energy from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("CAPE_mean_nr_E5",
       "ERA5 Rain-checked Mean CAPE",description,("tE5"),
       np.float64,CAPEE5units,fileout,CAPE_mean_nr_E5,
       f,float(-999))

#==================================================================
# Write SR18 information to file
#==================================================================

  if nl.addSR18E5=="True":

    if nl.addctarea=="True" and nl.addnorainchk=="True":
      description = "1000-850 hPa zonal shear from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("USR18_area_E5",
      "ERA5 1000-850 hPa zonal shear",description,
      ("tE5","yE5","xE5"),np.float64,SR18E5units,
       fileout,USR18_area_E5,f,float(-999))
      description = "1000-850 hPa meridional shear from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("VSR18_area_E5",
      "ERA5 1000-850 hPa meridional shear",description,
      ("tE5","yE5","xE5"),np.float64,SR18E5units,
       fileout,VSR18_area_E5,f,float(-999))

    if nl.addctmean=="True"and nl.addnorainchk=="True":
      description = "1000-850 hPa zonal shear from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("USR18_mean_E5",
       "ERA5 Mean 1000-850 hPa zonal shear",description,("tE5"),
       np.float64,SR18E5units,fileout,USR18_mean_E5,f,float(-999))
      description = "1000-850 hPa meridional shear from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("VSR18_mean_E5",
       "ERA5 Mean 1000-850 hPa Meridional Shear",description,
       ("tE5"),np.float64,SR18E5units,fileout,VSR18_mean_E5,
       f,float(-999))

    if nl.addctarea=="True"and nl.addrainchk=="True":
      USR18_area_nr_E5 = np.where(np.isnan(USR18_area_nr_E5),
       float(-999),USR18_area_nr_E5)
      description = "1000-850 hPa zonal shear from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("USR18_area_nr_E5",
       "ERA5 Rain-Checked 1000-850 Zonal Shear",description,
       ("tE5","yE5","xE5"),np.float64,SR18E5units,
       fileout,USR18_area_nr_E5,f,float(-999))
      description = "1000-850 hPa meridional shear from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("VSR18_area_nr_E5",
       "ERA5 Rain-Checked 1000-850 Meridional Shear",description,
       ("tE5","yE5","xE5"),np.float64,SR18E5units,
       fileout,VSR18_area_nr_E5,f,float(-999))

    if nl.addctmean=="True"and nl.addrainchk=="True":
      description = "1000-850 hPa zonal shear from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("USR18_mean_nr_E5",
       "ERA5 Rain-Checked Mean 1000-850 hPa Zonal Shear",
       description,("tE5"),np.float64,SR18E5units,fileout,
       USR18_mean_nr_E5,f,float(-999))
      description = "1000-850 hPa meridional shear from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("VSR18_mean_nr_E5",
       "ERA5 Rain-Checked Mean 1000-850 hPa Meridional Shear",
       description,("tE5"),np.float64,SR18E5units,fileout,
       VSR18_mean_nr_E5,f,float(-999))

#==================================================================
# Write SR82 information to file
#==================================================================

  if nl.addSR82E5=="True":

    if nl.addctarea=="True" and nl.addnorainchk=="True":
      description = "850-200 hPa zonal shear from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("USR82_area_E5",
      "ERA5 850-200 hPa zonal shear",description,
      ("tE5","yE5","xE5"),np.float64,SR82E5units,
       fileout,USR82_area_E5,f,float(-999))
      description = "850-200 hPa meridional shear from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("VSR82_area_E5",
      "ERA5 850-200 hPa meridional shear",description,
      ("tE5","yE5","xE5"),np.float64,SR82E5units,
       fileout,VSR82_area_E5,f,float(-999))

    if nl.addctmean=="True"and nl.addnorainchk=="True":
      description = "850-200 hPa zonal shear from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("USR82_mean_E5",
       "ERA5 Mean 850-200 hPa zonal shear",description,("tE5"),
       np.float64,SR82E5units,fileout,USR82_mean_E5,f,float(-999))
      description = "850-200 hPa meridional shear from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("VSR82_mean_E5",
       "ERA5 Mean 850-200 hPa Meridional Shear",description,
       ("tE5"),np.float64,SR82E5units,fileout,VSR82_mean_E5,
       f,float(-999))

    if nl.addctarea=="True"and nl.addrainchk=="True":
      USR82_area_nr_E5 = np.where(np.isnan(USR82_area_nr_E5),
       float(-999),USR82_area_nr_E5)
      description = "850-200 hPa zonal shear from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("USR82_area_nr_E5",
       "ERA5 Rain-Checked 850-200 Zonal Shear",description,
       ("tE5","yE5","xE5"),np.float64,SR82E5units,
       fileout,USR82_area_nr_E5,f,float(-999))
      description = "850-200 hPa meridional shear from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("VSR82_area_nr_E5",
       "ERA5 Rain-Checked 850-200 Meridional Shear",description,
       ("tE5","yE5","xE5"),np.float64,SR82E5units,
       fileout,VSR82_area_nr_E5,f,float(-999))

    if nl.addctmean=="True"and nl.addrainchk=="True":
      description = "850-200 hPa zonal shear from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("USR82_mean_nr_E5",
       "ERA5 Rain-Checked Mean 850-200 hPa Zonal Shear",
       description,("tE5"),np.float64,SR82E5units,fileout,
       USR82_mean_nr_E5,f,float(-999))
      description = "850-200 hPa meridional shear from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("VSR82_mean_nr_E5",
       "ERA5 Rain-Checked Mean 850-200 hPa Meridional Shear",
       description,("tE5"),np.float64,SR82E5units,fileout,
       VSR82_mean_nr_E5,f,float(-999))

#==================================================================
# Write SH18 information to file
#==================================================================

  if nl.addSH18E5=="True":

    if nl.addctarea=="True" and nl.addnorainchk=="True":
      description = "1000-850 hPa specific humidity from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("SH18_area_E5","ERA5 SH18",
       description,("tE5","yE5","xE5"),np.float64,SH18E5units,
       fileout,SH18_area_E5,f,float(-999))

    if nl.addctmean=="True"and nl.addnorainchk=="True":
      description = "1000-850 hPa specific humidity from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("SH18_mean_E5","ERA5 Mean SH18",
       description,("tE5"),np.float64,SH18E5units,fileout,
       SH18_mean_E5,f,float(-999))

    if nl.addctarea=="True"and nl.addrainchk=="True":
      SH18_area_nr_E5 = np.where(np.isnan(SH18_area_nr_E5),float(-999),SH18_area_nr_E5)
      description = "1000-850 hPa specific humidity from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("SH18_area_nr_E5","ERA5 Rain-Checked SH18",
       description,("tE5","yE5","xE5"),np.float64,SH18E5units,
       fileout,SH18_area_nr_E5,f,float(-999))

    if nl.addctmean=="True"and nl.addrainchk=="True":
      description = "1000-850 hPa specific humidity from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("SH18_mean_nr_E5",
       "ERA5 Rain-checked Mean SH18",description,("tE5"),
       np.float64,SH18E5units,fileout,SH18_mean_nr_E5,
       f,float(-999))

#==================================================================
# Write SH82 information to file
#==================================================================

  if nl.addSH82E5=="True":

    if nl.addctarea=="True" and nl.addnorainchk=="True":
      description = "850-200 hPa specific humidity from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("SH82_area_E5","ERA5 SH82",
       description,("tE5","yE5","xE5"),np.float64,SH82E5units,
       fileout,SH82_area_E5,f,float(-999))

    if nl.addctmean=="True"and nl.addnorainchk=="True":
      description = "850-200 hPa specific humidity from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("SH82_mean_E5","ERA5 Mean SH82",
       description,("tE5"),np.float64,SH82E5units,fileout,
       SH82_mean_E5,f,float(-999))

    if nl.addctarea=="True"and nl.addrainchk=="True":
      SH82_area_nr_E5 = np.where(np.isnan(SH82_area_nr_E5),float(-999),SH82_area_nr_E5)
      description = "850-200 hPa specific humidity from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("SH82_area_nr_E5","ERA5 Rain-Checked SH82",
       description,("tE5","yE5","xE5"),np.float64,SH82E5units,
       fileout,SH82_area_nr_E5,f,float(-999))

    if nl.addctmean=="True"and nl.addrainchk=="True":
      description = "850-200 hPa specific humidity from the ERA5 dataset for a 10x10 degree area centered on the precipitation system centroid. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
      mfns.write_var("SH82_mean_nr_E5",
       "ERA5 Rain-checked Mean SH82",description,("tE5"),
       np.float64,SH82E5units,fileout,SH82_mean_nr_E5,
       f,float(-999))

#==================================================================
# Close current file
#==================================================================

  fileout.close()

#==================================================================
# End processing of current PF
#==================================================================
