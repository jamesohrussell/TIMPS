#==================================================================
# Add ERA5 variables to IMERG PF netcdf files. Requires 
#  driver_addERA5vars.py be in the same directory.
#
# James Russell 2019
#==================================================================

# Import python libraries (do not change)
import numpy as np
from netCDF4 import Dataset
import glob
import datetime as dt
from joblib import Parallel, delayed
import time as tm
import os
import sys

# Import namelist
from namelist_TIMPS import general as gnl
from namelist_TIMPS import addERA5vars as anl

# Import custom libraries
sys.path.insert(0,gnl.fnsdir)
import time_functions as tfns
import misc_functions as mfns
import ERA5_functions as E5fns

import warnings
warnings.filterwarnings("ignore")

#==================================================================
# Initialize timer
#==================================================================

startnow = tm.time()

#==================================================================
# Function to pull and process data
#==================================================================

def create_var(varname,varfname,files,timestr,anl,TCRW,loni,lati):

  # Find file and time indices
  fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       files[varname],timestr)

  # Get an area
  if anl.addfull:
    var = E5fns.get_E5_ss_2D_var(
     fh,varfname,timi,loni,lati,times,ctime)
    var_nr = np.where(TCRW>0.001,np.nan,var)
    varnan = var_nr.flatten()
    if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
      var_mean_nr = np.nanmean(var_nr)
    else:
      var_mean_nr = np.nan
  if anl.addanom:
    varanom = E5fns.get_E5_ss_2D_var(
     fh,f'{varfname}_anom',timi,loni,lati,times,ctime)
    varanom_nr = np.where(TCRW>0.001,np.nan,varanom)
    varnan = varanom_nr.flatten()
    if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
      varanom_mean_nr = np.nanmean(varanom_nr)
    else:
      varanom_mean_nr = np.nan
  if anl.addanomn:
    varanomn = E5fns.get_E5_ss_2D_var(
     fh,f'{varfname}_anomn',timi,loni,lati,times,ctime)
    varanomn_nr = np.where(TCRW>0.001,np.nan,varanomn)
    varnan = varanomn_nr.flatten()
    if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
      varanomn_mean_nr = np.nanmean(varanomn_nr)
    else:
      varanomn_mean_nr = np.nan

  # Close file
  for fii in fh:
    fii.close()

  return(var_mean_nr,varanom_mean_nr,varanomn_mean_nr)

#==================================================================
# Function to write data
#==================================================================

def write_data(lname,sname,fileg,var,vara,varan,units):

    description = f"{lname} from the ERA5 dataset averaged over a 5x5 degree area, centered on the precipitation system centroid, and areas with substantial rain water removed. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    if anl.addfull:
      mfns.write_var(f"{sname}_mean_nr",
       f"ERA5 Mean {lname} (Rain cells removed)",
       description,("tE5"),np.float64,units,fileg,var)
    if anl.addanom:
      mfns.write_var(f"{sname}_anom_mean_nr",
       f"ERA5 Mean {lname} Anomaly (Rain cells removed)",
       description,("tE5"),np.float64,units,fileg,vara)
    if anl.addanomn:
      mfns.write_var(f"{sname}_anomn_mean_nr",
       f"ERA5 Mean {lname} Normalized Anomaly (Rain cells removed)",
       description,("tE5"),np.float64,"",fileg,varan)

#==================================================================
# Generate list of TIMPS files to process
#==================================================================

# Reads directory and file names
if anl.ssdat:

  print("Subsetting by date")

  # Generate a list of filenames with dates to search for
  print("Generating filenames to search for")
  start = dt.datetime.strptime(anl.date1,"%Y%m%d")
  end = dt.datetime.strptime(anl.date2,"%Y%m%d")
  datearr = (start + dt.timedelta(days=x) \
   for x in range(0,(end-start).days))
  filen = []
  for dateobj in datearr:
    filen.extend(f'{gnl.datadirTIPS}{dateobj.strftime("%m")}/\
{gnl.fileidTIPS}*_{dateobj.strftime("%Y%m%d")}*.nc')

  # Reads directory and filenames
  print("Getting directory and filenames")
  filenamesrun = []
  for f in filen:
    filenamesrun.extend(sorted(glob.glob(f)))

elif anl.ssobs:

  print("Subsetting by object ids")

  print("Generating filenames to search for")
  obids = [i for i in range(int(anl.obid1),int(anl.obid2))]
  filen = []
  for j in range(1,13):
    for i in range(len(obids)): 
      filen.extend(f'{gnl.datadirTIPS}{str(j).zfill(2)}/\
{gnl.fileidTIPS}{str(obids[i])}*')

  # Reads directory and filenames
  print("Getting directory and filenames")
  filenamesrun = []
  for f in filen:
    print(f'Adding {str(f)}')
    filenamesrun.extend(sorted(glob.glob(f)))

elif anl.sslst:

  print("Subsetting by files in csv")

  print("Generating list of files")
  filenamesrun = list(np.genfromtxt(anl.lstfn,
   delimiter=',',dtype=str))

else:

  print("No subsetting")

  # Find all files in directory
  print("Generating list of files")
  filenamesrun = []
  for i in range(1,13):
    filenamesrun.extend(sorted(glob.glob(
     f'{gnl.datadirTIPS}{str(i).zfill(2)}/{gnl.fileidTIPS}*')))

#==================================================================
# Define function
#==================================================================

def driver_addvars(fn):
  """
  Driver to calculate various variables for a single PF.
  fn corresponds to the line number in filenames_av.txt.
  """
  warnings.filterwarnings("ignore")

#==================================================================
# Begin script
#==================================================================

  # Read in filename for PF
  f = filenamesrun[fn]
  print(f'Working with file {str(f)}')

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
  timestrs = [f'{str(it)[0:4]}-{str(it)[4:6]}-{str(it)[6:8]} \
{str(it)[8:10]}:{str(it)[10:12]}:00' for it in datadtim]

  # Define datetime objects for first and last times
  fto = dt.datetime(int(timestrs[0][0:4]),int(timestrs[0][5:7]),
   int(timestrs[0][8:10]),hour=int(timestrs[0][11:13]),
   minute=int(timestrs[0][14:16]),second=int(timestrs[0][17:19]))
  lto = dt.datetime(int(timestrs[-1][0:4]),int(timestrs[-1][5:7]),
   int(timestrs[-1][8:10]),hour=int(timestrs[-1][11:13]),
   minute=int(timestrs[-1][14:16]),second=int(timestrs[-1][17:19]))
  
  # Find times before and after and add to time list
  timestrs = [str(fto-dt.timedelta(hours=int(it))) 
              for it in anl.hoursbefore] + timestrs + \
             [str(lto+dt.timedelta(hours=int(it))) 
              for it in anl.hoursafter]
  tE5 = [tfns.time_since(i,timeunits) for i in timestrs]

  # Find files
  files["TCRW"] = E5fns.get_E5_ss_files(
    anl.fileTCRWid,timestrs[0],timestrs[-1])

  # Moisture variables
  if anl.addTCWV: files["TCWV"] = E5fns.get_E5_ss_files(
    anl.fileTCWVid,timestrs[0],timestrs[-1])
  if anl.addSH18: files["SH18"] = E5fns.get_E5_ss_files(
    anl.fileSH18id,timestrs[0],timestrs[-1])
  if anl.addSH84: files["SH84"] = E5fns.get_E5_ss_files(
    anl.fileSH84id,timestrs[0],timestrs[-1])
  if anl.addMF18: files["MF18"] = E5fns.get_E5_ss_files(
    anl.fileMF18id,timestrs[0],timestrs[-1])
  if anl.addMF84: files["MF84"] = E5fns.get_E5_ss_files(
    anl.fileMF84id,timestrs[0],timestrs[-1])
  if anl.addMA18: files["MA18"] = E5fns.get_E5_ss_files(
    anl.fileMA18id,timestrs[0],timestrs[-1])
  if anl.addMA84: files["MA84"] = E5fns.get_E5_ss_files(
    anl.fileMA84id,timestrs[0],timestrs[-1])
  if anl.addMC18: files["MC18"] = E5fns.get_E5_ss_files(
    anl.fileMC18id,timestrs[0],timestrs[-1])
  if anl.addMC84: files["MC84"] = E5fns.get_E5_ss_files(
    anl.fileMC84id,timestrs[0],timestrs[-1])

  # Kinematic variables
  if anl.addSR18: files["SR18"] = E5fns.get_E5_ss_files(
    anl.fileSR18id,timestrs[0],timestrs[-1])
  if anl.addSR17: files["SR17"] = E5fns.get_E5_ss_files(
    anl.fileSR18id,timestrs[0],timestrs[-1])
  if anl.addSR84: files["SR84"] = E5fns.get_E5_ss_files(
    anl.fileSR84id,timestrs[0],timestrs[-1])
  if anl.addSR65: files["SR65"] = E5fns.get_E5_ss_files(
    anl.fileSR65id,timestrs[0],timestrs[-1])
  if anl.addSR14: files["SR14"] = E5fns.get_E5_ss_files(
    anl.fileSR14id,timestrs[0],timestrs[-1])
  if anl.addCV18: files["CV18"] = E5fns.get_E5_ss_files(
    anl.fileCV18id,timestrs[0],timestrs[-1])
  if anl.addCV31: files["CV31"] = E5fns.get_E5_ss_files(
    anl.fileCV31id,timestrs[0],timestrs[-1])
  if anl.addVE18: files["VE18"] = E5fns.get_E5_ss_files(
    anl.fileVE18id,timestrs[0],timestrs[-1])
  if anl.addVE84: files["VE84"] = E5fns.get_E5_ss_files(
    anl.fileVE84id,timestrs[0],timestrs[-1])

#==================================================================
# Get coordinates and their indices
#==================================================================

  # Get corresponding lists of central latitude and longitude
  dataclat = [round(i,2) for i in 
             [dataclat[0]]*len(anl.hoursbefore) + \
             [i for i in dataclat] + \
             [dataclat[-1]]*len(anl.hoursafter)]
  dataclon = [round(i,2) for i in 
             [dataclon[0]]*len(anl.hoursbefore) + \
             [i for i in dataclon] + \
             [dataclon[-1]]*len(anl.hoursafter)]

#==================================================================
# Begin loop over times
#==================================================================

  c = 0
  for k in timestrs:

#==================================================================
# Get coordinates
#==================================================================

    # Find coordinates and indices
    keys = [str(k) for k in files.keys()]
    loni,lati = E5fns.get_E5_ss_coords(
        Dataset(files[keys[0]][0]),dataclon[c],dataclat[c],anl.hda)

#==================================================================
# Preallocate arrays
#==================================================================

    if c==0:

      # Moisture

      if anl.addTCWV:

        if anl.addfull:
          TCWV_mean_nr = [np.nan]*len(timestrs)
        if anl.addanom:
          TCWVanom_mean_nr = [np.nan]*len(timestrs)
        if anl.addanomn:
          TCWVanomn_mean_nr = [np.nan]*len(timestrs)

      if anl.addSH18:

        if anl.addfull:
          SH18_mean_nr = [np.nan]*len(timestrs)
        if anl.addanom:
          SH18anom_mean_nr = [np.nan]*len(timestrs)
        if anl.addanomn:
          SH18anomn_mean_nr = [np.nan]*len(timestrs)

      if anl.addSH84:

        if anl.addfull:
          SH84_mean_nr = [np.nan]*len(timestrs)
        if anl.addanom:
          SH84anom_mean_nr = [np.nan]*len(timestrs)
        if anl.addanomn:
          SH84anomn_mean_nr = [np.nan]*len(timestrs)

      if anl.addMF18:

        if anl.addfull:
          MF18_mean_nr = [np.nan]*len(timestrs)
        if anl.addanom:
          MF18anom_mean_nr = [np.nan]*len(timestrs)
        if anl.addanomn:
          MF18anomn_mean_nr = [np.nan]*len(timestrs)

      if anl.addMF84:

        if anl.addfull:
          MF84_mean_nr = [np.nan]*len(timestrs)
        if anl.addanom:
          MF84anom_mean_nr = [np.nan]*len(timestrs)
        if anl.addanomn:
          MF84anomn_mean_nr = [np.nan]*len(timestrs)

      if anl.addMA18:

        if anl.addfull:
          MA18_mean_nr = [np.nan]*len(timestrs)
        if anl.addanom:
          MA18anom_mean_nr = [np.nan]*len(timestrs)
        if anl.addanomn:
          MA18anomn_mean_nr = [np.nan]*len(timestrs)

      if anl.addMA84:

        if anl.addfull:
          MA84_mean_nr = [np.nan]*len(timestrs)
        if anl.addanom:
          MA84anom_mean_nr = [np.nan]*len(timestrs)
        if anl.addanomn:
          MA84anomn_mean_nr = [np.nan]*len(timestrs)

      if anl.addMC18:

        if anl.addfull:
          MC18_mean_nr = [np.nan]*len(timestrs)
        if anl.addanom:
          MC18anom_mean_nr = [np.nan]*len(timestrs)
        if anl.addanomn:
          MC18anomn_mean_nr = [np.nan]*len(timestrs)

      if anl.addMC84:

        if anl.addfull:
          MC84_mean_nr = [np.nan]*len(timestrs)
        if anl.addanom:
          MC84anom_mean_nr = [np.nan]*len(timestrs)
        if anl.addanomn:
          MC84anomn_mean_nr = [np.nan]*len(timestrs)

      # Shear

      if anl.addSR18:

        if anl.addfull:
          SR18_mean_nr = [np.nan]*len(timestrs)
        if anl.addanom:
          SR18anom_mean_nr = [np.nan]*len(timestrs)
        if anl.addanomn:
          SR18anomn_mean_nr = [np.nan]*len(timestrs)

      if anl.addSR17:

        if anl.addfull:
          SR17_mean_nr = [np.nan]*len(timestrs)
        if anl.addanom:
          SR17anom_mean_nr = [np.nan]*len(timestrs)
        if anl.addanomn:
          SR17anomn_mean_nr = [np.nan]*len(timestrs)

      if anl.addSR84:

        if anl.addfull:
          SR84_mean_nr = [np.nan]*len(timestrs)
        if anl.addanom:
          SR84anom_mean_nr = [np.nan]*len(timestrs)
        if anl.addanomn:
          SR84anomn_mean_nr = [np.nan]*len(timestrs)

      if anl.addSR65:

        if anl.addfull:
          SR65_mean_nr = [np.nan]*len(timestrs)
        if anl.addanom:
          SR65anom_mean_nr = [np.nan]*len(timestrs)
        if anl.addanomn:
          SR65anomn_mean_nr = [np.nan]*len(timestrs)

      if anl.addSR14:

        if anl.addfull:
          SR14_mean_nr = [np.nan]*len(timestrs)
        if anl.addanom:
          SR14anom_mean_nr = [np.nan]*len(timestrs)
        if anl.addanomn:
          SR14anomn_mean_nr = [np.nan]*len(timestrs)

      # Convergence/divergence

      if anl.addCV18:

        if anl.addfull:
          CV18_mean_nr = [np.nan]*len(timestrs)
        if anl.addanom:
          CV18anom_mean_nr = [np.nan]*len(timestrs)
        if anl.addanomn:
          CV18anomn_mean_nr = [np.nan]*len(timestrs)

      if anl.addCV31:

        if anl.addfull:
          CV31_mean_nr = [np.nan]*len(timestrs)
        if anl.addanom:
          CV31anom_mean_nr = [np.nan]*len(timestrs)
        if anl.addanomn:
          CV31anomn_mean_nr = [np.nan]*len(timestrs)

      # Vertical motion
      if anl.addVE18:

        if anl.addfull:
          VE18_mean_nr = [np.nan]*len(timestrs)
        if anl.addanom:
          VE18anom_mean_nr = [np.nan]*len(timestrs)
        if anl.addanomn:
          VE18anomn_mean_nr = [np.nan]*len(timestrs)

      if anl.addVE84:

        if anl.addfull:
          VE84_mean_nr = [np.nan]*len(timestrs)
        if anl.addanom:
          VE84anom_mean_nr = [np.nan]*len(timestrs)
        if anl.addanomn:
          VE84anomn_mean_nr = [np.nan]*len(timestrs)

#==================================================================
# Get rain check data
#==================================================================

    if anl.addrainchk:
      # Find file and time indices
      fhR,timiR,timesR,ctimeR = E5fns.get_E5_ss_2D_fiti(
       files["TCRW"],timestrs[c])

      TCRW = E5fns.get_E5_ss_2D_var(
       fhR,"TCRW",timiR,loni,lati,timesR,ctimeR)

      # Close file
      for fii in fhR:
        fii.close()

#==================================================================
# Get data
#==================================================================

    # TCWV
    if anl.addTCWV:
      TCWV_mean_nr[c],TCWVanom_mean_nr[c],TCWVanomn_mean_nr[c] = \
       create_var("TCWV","TCWV",files,timestrs[c],anl,TCRW,
       loni,lati)

    # Specific humidity
    if anl.addSH18:
      SH18_mean_nr[c],SH18anom_mean_nr[c],SH18anomn_mean_nr[c] = \
       create_var("SH18","Q",files,timestrs[c],anl,TCRW,
       loni,lati)

    if anl.addSH84:
      SH84_mean_nr[c],SH84anom_mean_nr[c],SH84anomn_mean_nr[c] = \
       create_var("SH84","Q",files,timestrs[c],anl,TCRW,
       loni,lati)

    # Moisture flux
    if anl.addMF18:
      MF18_mean_nr[c],MF18anom_mean_nr[c],MF18anomn_mean_nr[c] = \
       create_var("MF18","MFC",files,timestrs[c],anl,TCRW,
       loni,lati)

    if anl.addMF84:
      MF84_mean_nr[c],MF84anom_mean_nr[c],MF84anomn_mean_nr[c] = \
       create_var("MF84","MFC",files,timestrs[c],anl,TCRW,
       loni,lati)

    # Moisture advection
    if anl.addMA18:
      MA18_mean_nr[c],MA18anom_mean_nr[c],MA18anomn_mean_nr[c] = \
       create_var("MA18","MADV",files,timestrs[c],anl,TCRW,
       loni,lati)

    if anl.addMA84:
      MA84_mean_nr[c],MA84anom_mean_nr[c],MA84anomn_mean_nr[c] = \
       create_var("MA84","MADV",files,timestrs[c],anl,TCRW,
       loni,lati)

    # Moisture convergence
    if anl.addMC18:
      MC18_mean_nr[c],MC18anom_mean_nr[c],MC18anomn_mean_nr[c] = \
       create_var("MC18","MCONV",files,timestrs[c],anl,TCRW,
       loni,lati)

    if anl.addMC84:
      MC84_mean_nr[c],MC84anom_mean_nr[c],MC84anomn_mean_nr[c] = \
       create_var("MC84","MCONV",files,timestrs[c],anl,TCRW,
       loni,lati)

    # Shear
    if anl.addSR18:
      SR18_mean_nr[c],SR18anom_mean_nr[c],SR18anomn_mean_nr[c] = \
       create_var("SR18","MSHR",files,timestrs[c],anl,TCRW,
       loni,lati)

    if anl.addSR17:
      SR17_mean_nr[c],SR17anom_mean_nr[c],SR17anomn_mean_nr[c] = \
       create_var("SR17","MSHR",files,timestrs[c],anl,TCRW,
       loni,lati)

    if anl.addSR84:
      SR84_mean_nr[c],SR84anom_mean_nr[c],SR84anomn_mean_nr[c] = \
       create_var("SR84","MSHR",files,timestrs[c],anl,TCRW,
       loni,lati)

    if anl.addSR65:
      SR65_mean_nr[c],SR65anom_mean_nr[c],SR65anomn_mean_nr[c] = \
       create_var("SR65","MSHR",files,timestrs[c],anl,TCRW,
       loni,lati)

    if anl.addSR14:
      SR14_mean_nr[c],SR14anom_mean_nr[c],SR14anomn_mean_nr[c] = \
       create_var("SR14","MSHR",files,timestrs[c],anl,TCRW,
       loni,lati)

    # Convergence
    if anl.addCV18:
      CV18_mean_nr[c],CV18anom_mean_nr[c],CV18anomn_mean_nr[c] = \
       create_var("CV18","CONV",files,timestrs[c],anl,TCRW,
       loni,lati)

    if anl.addCV31:
      CV31_mean_nr[c],CV31anom_mean_nr[c],CV31anomn_mean_nr[c] = \
       create_var("CV31","CONV",files,timestrs[c],anl,TCRW,
       loni,lati)

    # Vertical motion
    if anl.addVE18:
      VE18_mean_nr[c],VE18anom_mean_nr[c],VE18anomn_mean_nr[c] = \
       create_var("VE18","W",files,timestrs[c],anl,TCRW,
       loni,lati)

    if anl.addVE84:
      VE84_mean_nr[c],VE84anom_mean_nr[c],VE84anomn_mean_nr[c] = \
       create_var("VE84","W",files,timestrs[c],anl,TCRW,
       loni,lati)

#==================================================================
# End loops over objects and times
#==================================================================
      
    # Advance counter
    c = c + 1

  fd.close()

#==================================================================
# Open file to write data and create ERA5 group
#==================================================================

  # Open file
  fileout = Dataset(f,'a')

  # Create group
  try:
    E5group  = fileout.createGroup("ERA5")
  except:
    E5group = fileout.group["ERA5"]

#==================================================================
# Write coordinate data for environment information to file
#==================================================================

  try: tE51 = E5group.createDimension('tE5',len(tE5))
  except: print("tE5 already defined")

  description = "Corresponds to ERA5 data. Different from time dimension since times are added before and after the TIPS exists."
  mfns.write_var("tE5","ERA5 time",description,("tE5"),
    np.float64,timeunits,E5group,tE5)

#==================================================================
# Write TCWV information to file
#==================================================================

  # TCWV
  if anl.addTCWV:
    write_data("Total column water vapor","TCWV",E5group,
     TCWV_mean_nr,TCWVanom_mean_nr,TCWVanomn_mean_nr,"kg m**-2")

  # Specific humidity
  if anl.addSH18:
    write_data("Specific humidity (1000-850hPa average)","SH18",
     E5group,SH18_mean_nr,SH18anom_mean_nr,SH18anomn_mean_nr,
     "g kg**-1")

  if anl.addSH84:
    write_data("Specific humidity (800-400hPa average)","SH84",
     E5group,SH84_mean_nr,SH84anom_mean_nr,SH84anomn_mean_nr,
     "g kg**-1")

  # Moisture flux
  if anl.addMF18:
    write_data("Moisture flux convergence (1000-850hPa average)",
     "MFC18",E5group,MF18_mean_nr,MF18anom_mean_nr,
     MF18anomn_mean_nr,"g kg**-1 hr**-1")

  if anl.addMF84:
    write_data("Moisture flux convergence (800-400hPa average)",
     "MFC84",E5group,MF84_mean_nr,MF84anom_mean_nr,
     MF84anomn_mean_nr,"g kg**-1 hr**-1")

  # Moisture advection
  if anl.addMA18:
    write_data("Moisture advection (1000-850hPa average)",
     "MADV18",E5group,MA18_mean_nr,MA18anom_mean_nr,
     MA18anomn_mean_nr,"g kg**-1 hr**-1")

  if anl.addMA84:
    write_data("Moisture advection (800-400hPa average)",
     "MADV84",E5group,MA84_mean_nr,MA84anom_mean_nr,
     MA84anomn_mean_nr,"g kg**-1 hr**-1")

  # Moisture convergence
  if anl.addMC18:
    write_data("Moisture convergence (1000-850hPa average)",
     "MCONV18",E5group,MC18_mean_nr,MC18anom_mean_nr,
     MC18anomn_mean_nr,"g kg**-1 hr**-1")

  if anl.addMC84:
    write_data("Moisture convergence (800-400hPa average)",
     "MCONV84",E5group,MC84_mean_nr,MC84anom_mean_nr,
     MC84anomn_mean_nr,"g kg**-1 hr**-1")

  # Shear
  if anl.addSR18:
    write_data("Shear magnitude (1000-850hPa average, ~0-1.5km)",
     "MSHR18",E5group,SR18_mean_nr,SR18anom_mean_nr,
     SR18anomn_mean_nr,"m s**-1")

  if anl.addSR17:
    write_data("Shear magnitude (1000-700hPa average, ~0-3km)",
     "MSHR17",E5group,SR17_mean_nr,SR17anom_mean_nr,
     SR17anomn_mean_nr,"m s**-1")

  if anl.addSR84:
    write_data("Shear magnitude (800-400hPa average, ~2-8km)",
     "MSHR84",E5group,SR84_mean_nr,SR84anom_mean_nr,
     SR84anomn_mean_nr,"m s**-1")

  if anl.addSR65:
    write_data("Shear magnitude (650-500hPa average, ~4-6km)",
     "MSHR65",E5group,SR65_mean_nr,SR65anom_mean_nr,
     SR65anomn_mean_nr,"m s**-1")

  if anl.addSR14:
    write_data("Shear magnitude (1000-400hPa average, ~0-8km)",
     "MSHR14",E5group,SR14_mean_nr,SR14anom_mean_nr,
     SR14anomn_mean_nr,"m s**-1")

  # Convergence
  if anl.addCV18:
    write_data("Convergence (1000-850hPa average)",
     "CONV18",E5group,CV18_mean_nr,CV18anom_mean_nr,
     CV18anomn_mean_nr,"hr**-1")

  if anl.addCV31:
    write_data("Convergence (300-100hPa average)",
     "CONV31",E5group,CV31_mean_nr,CV31anom_mean_nr,
     CV31anomn_mean_nr,"hr**-1")

  # Vertical motion
  if anl.addVE18:
    write_data("Vertical motion (1000-850hPa average)",
     "VERT18",E5group,VE18_mean_nr,VE18anom_mean_nr,
     VE18anomn_mean_nr,"Pa s**-1")

  if anl.addVE84:
    write_data("Vertical motion (800-400hPa average)",
     "VERT84",E5group,VE84_mean_nr,VE84anom_mean_nr,
     VE84anomn_mean_nr,"Pa s**-1")

#==================================================================
# Close current file
#==================================================================

  fileout.close()

#==================================================================
# End processing of current PF
#==================================================================

#==================================================================
# Loop over IPF files in parrallel
#==================================================================

if anl.serialorparallel==1:
  print("Begin serial loop over TIMPS")
  for fn in range(len(filenamesrun)):
    driver_addvars(fn)

# Parrallel loop over PFs
if anl.serialorparallel==2:
  print("Begin parrallel loop over TIMPS")
  Parallel(n_jobs=anl.njobs)(delayed(driver_addvars)(fn) \
    for fn in range(len(filenamesrun)))

#==================================================================
# Final clean up tasks
#==================================================================

# Remove namelist and filename files
print("Cleaning up")
os.system("rm -rf __pycache__")

# End timer
endnow = tm.time()
print(f'Program took {str(endnow - startnow)}s')

#==================================================================
# End code
#==================================================================
