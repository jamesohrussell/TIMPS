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
    anl.fileTCRWE5id,timestrs[0],timestrs[-1])

  # Moisture variables
  if anl.addTCWVE5: files["TCWV"] = E5fns.get_E5_ss_files(
    anl.fileTCWVE5id,timestrs[0],timestrs[-1])
  if anl.addSH18E5: files["SH18"] = E5fns.get_E5_ss_files(
    anl.fileSH18E5id,timestrs[0],timestrs[-1])
  if anl.addSH84E5: files["SH84"] = E5fns.get_E5_ss_files(
    anl.fileSH84E5id,timestrs[0],timestrs[-1])

  # Shear variables
  if anl.addSR18E5: files["SR18"] = E5fns.get_E5_ss_files(
    anl.fileSR18E5id,timestrs[0],timestrs[-1])
  if anl.addSR84E5: files["SR84"] = E5fns.get_E5_ss_files(
    anl.fileSR84E5id,timestrs[0],timestrs[-1])
  if anl.addSR14E5: files["SR14"] = E5fns.get_E5_ss_files(
    anl.fileSR14E5id,timestrs[0],timestrs[-1])

  # Convergence/divergence variables
  if anl.addCV18E5: files["CV18"] = E5fns.get_E5_ss_files(
    anl.fileCV18E5id,timestrs[0],timestrs[-1])
  if anl.addCV31E5: files["CV31"] = E5fns.get_E5_ss_files(
    anl.fileCV31E5id,timestrs[0],timestrs[-1])

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

      if anl.addTCWVE5:

        if anl.addfull:
          TCWV_mean_nr_E5 = [np.nan]*len(timestrs)
        if anl.addanom:
          TCWVanom_mean_nr_E5 = [np.nan]*len(timestrs)
        if anl.addanomn:
          TCWVanomn_mean_nr_E5 = [np.nan]*len(timestrs)

      if anl.addSH18E5:

        if anl.addfull:
          SH18_mean_nr_E5 = [np.nan]*len(timestrs)
        if anl.addanom:
          SH18anom_mean_nr_E5 = [np.nan]*len(timestrs)
        if anl.addanomn:
          SH18anomn_mean_nr_E5 = [np.nan]*len(timestrs)

      if anl.addSH84E5:

        if anl.addfull:
          SH84_mean_nr_E5 = [np.nan]*len(timestrs)
        if anl.addanom:
          SH84anom_mean_nr_E5 = [np.nan]*len(timestrs)
        if anl.addanomn:
          SH84anomn_mean_nr_E5 = [np.nan]*len(timestrs)

      # Shear

      if anl.addSR18E5:

        if anl.addfull:
          SR18_mean_nr_E5 = [np.nan]*len(timestrs)
        if anl.addanom:
          SR18anom_mean_nr_E5 = [np.nan]*len(timestrs)
        if anl.addanomn:
          SR18anomn_mean_nr_E5 = [np.nan]*len(timestrs)

      if anl.addSR84E5:

        if anl.addfull:
          SR84_mean_nr_E5 = [np.nan]*len(timestrs)
        if anl.addanom:
          SR84anom_mean_nr_E5 = [np.nan]*len(timestrs)
        if anl.addanomn:
          SR84anomn_mean_nr_E5 = [np.nan]*len(timestrs)

      if anl.addSR14E5:

        if anl.addfull:
          SR14_mean_nr_E5 = [np.nan]*len(timestrs)
        if anl.addanom:
          SR14anom_mean_nr_E5 = [np.nan]*len(timestrs)
        if anl.addanomn:
          SR14anomn_mean_nr_E5 = [np.nan]*len(timestrs)

      # Convergence/divergence

      if anl.addCV18E5:

        if anl.addfull:
          CV18_mean_nr_E5 = [np.nan]*len(timestrs)
        if anl.addanom:
          CV18anom_mean_nr_E5 = [np.nan]*len(timestrs)
        if anl.addanomn:
          CV18anomn_mean_nr_E5 = [np.nan]*len(timestrs)

      if anl.addCV31E5:

        if anl.addfull:
          CV31_mean_nr_E5 = [np.nan]*len(timestrs)
        if anl.addanom:
          CV31anom_mean_nr_E5 = [np.nan]*len(timestrs)
        if anl.addanomn:
          CV31anomn_mean_nr_E5 = [np.nan]*len(timestrs)

#==================================================================
# Get rain check data
#==================================================================

    if anl.addrainchk:
      # Find file and time indices
      fhR,timiR,timesR,ctimeR = E5fns.get_E5_ss_2D_fiti(
       files["TCRW"],timestrs[c])

      TCRWE5 = E5fns.get_E5_ss_2D_var(
       fhR,"TCRW",timiR,loni,lati,timesR,ctimeR)

      # Close file
      for fii in fhR:
        fii.close()

#==================================================================
# Get TCWV data
#==================================================================

    if anl.addTCWVE5:

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       files["TCWV"],timestrs[c])

      # Get units
      TCWVE5units = fh[0].variables["TCWV"].units

      # Get an area
      if anl.addfull:
        TCWVE5 = E5fns.get_E5_ss_2D_var(
               fh,"TCWV",timi,loni,lati,times,ctime)
        TCWVE5_nr = np.where(TCRWE5>0.001,np.nan,TCWVE5)
        varnan = TCWVE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          TCWV_mean_nr_E5[c] = np.nanmean(TCWVE5_nr)

      if anl.addanom:
        TCWVanomE5 = E5fns.get_E5_ss_2D_var(
               fh,"TCWV_anom",timi,loni,lati,times,ctime)
        TCWVanomE5_nr = np.where(TCRWE5>0.001,np.nan,TCWVanomE5)
        varnan = TCWVanomE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          TCWVanom_mean_nr_E5[c] = np.nanmean(TCWVanomE5_nr)
      if anl.addanomn:
        TCWVanomnE5 = E5fns.get_E5_ss_2D_var(
               fh,"TCWV_anomn",timi,loni,lati,times,ctime)
        TCWVanomnE5_nr = np.where(TCRWE5>0.001,np.nan,TCWVanomnE5)
        varnan = TCWVanomnE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          TCWVanomn_mean_nr_E5[c] = np.nanmean(TCWVanomnE5_nr)         
          
      # Close file
      for fii in fh:
        fii.close()

#==================================================================
# Get SH18 data
#==================================================================

    if anl.addSH18E5:

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       files["SH18"],timestrs[c])

      # Get units
      SH18E5units = fh[0].variables["Q"].units

      # Get an area
      if anl.addfull:
        SH18E5 = E5fns.get_E5_ss_2D_var(
               fh,"Q",timi,loni,lati,times,ctime)
        SH18E5_nr = np.where(TCRWE5>0.001,np.nan,SH18E5)
        varnan = SH18E5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          SH18_mean_nr_E5[c] = np.nanmean(SH18E5_nr)
      if anl.addanom:
        SH18anomE5 = E5fns.get_E5_ss_2D_var(
               fh,"Q_anom",timi,loni,lati,times,ctime)
        SH18anomE5_nr = np.where(TCRWE5>0.001,np.nan,SH18anomE5)
        varnan = SH18anomE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          SH18anom_mean_nr_E5[c] = np.nanmean(SH18anomE5_nr)
      if anl.addanomn:
        SH18anomnE5 = E5fns.get_E5_ss_2D_var(
               fh,"Q_anomn",timi,loni,lati,times,ctime)
        SH18anomnE5_nr = np.where(TCRWE5>0.001,np.nan,SH18anomnE5)
        varnan = SH18anomnE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          SH18anomn_mean_nr_E5[c] = np.nanmean(SH18anomnE5_nr)

      # Close file
      for fii in fh:
        fii.close()
     
#==================================================================
# Get SH84 data
#==================================================================

    if anl.addSH84E5:

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       files["SH84"],timestrs[c])

      # Get units
      SH84E5units = fh[0].variables["Q"].units

      # Get an area
      if anl.addfull:
        SH84E5 = E5fns.get_E5_ss_2D_var(
               fh,"Q",timi,loni,lati,times,ctime)
        SH84E5_nr = np.where(TCRWE5>0.001,np.nan,SH84E5)
        varnan = SH84E5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          SH84_mean_nr_E5[c] = np.nanmean(SH84E5_nr)
      if anl.addanom:
        SH84anomE5 = E5fns.get_E5_ss_2D_var(
               fh,"Q_anom",timi,loni,lati,times,ctime)
        SH84anomE5_nr = np.where(TCRWE5>0.001,np.nan,SH84anomE5)
        varnan = SH84anomE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          SH84anom_mean_nr_E5[c] = np.nanmean(SH84anomE5_nr)
      if anl.addanomn:
        SH84anomnE5 = E5fns.get_E5_ss_2D_var(
               fh,"Q_anomn",timi,loni,lati,times,ctime)
        SH84anomnE5_nr = np.where(TCRWE5>0.001,np.nan,SH84anomnE5)
        varnan = SH84anomnE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          SH84anomn_mean_nr_E5[c] = np.nanmean(SH84anomnE5_nr)

      # Close file
      for fii in fh:
        fii.close()

#==================================================================
# Get SR18 data
#==================================================================

    if anl.addSR18E5:

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       files["SR18"],timestrs[c])

      # Get units
      SR18E5units = fh[0].variables["MSHR"].units

      # Get an area
      if anl.addfull:
        SR18E5 = E5fns.get_E5_ss_2D_var(
               fh,"MSHR",timi,loni,lati,times,ctime)
        SR18E5_nr = np.where(TCRWE5>0.001,np.nan,SR18E5)
        varnan = SR18E5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          SR18_mean_nr_E5[c] = np.nanmean(SR18E5_nr)
      if anl.addanom:
        SR18anomE5 = E5fns.get_E5_ss_2D_var(
               fh,"MSHR_anom",timi,loni,lati,times,ctime)
        SR18anomE5_nr = np.where(TCRWE5>0.001,np.nan,SR18anomE5)
        varnan = SR18anomE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          SR18anom_mean_nr_E5[c] = np.nanmean(SR18anomE5_nr)
      if anl.addanomn:
        SR18anomnE5 = E5fns.get_E5_ss_2D_var(
               fh,"MSHR_anomn",timi,loni,lati,times,ctime)
        SR18anomnE5_nr = np.where(TCRWE5>0.001,np.nan,SR18anomnE5)
        varnan = SR18anomnE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          SR18anomn_mean_nr_E5[c] = np.nanmean(SR18anomnE5_nr)

      # Close file
      for fii in fh:
        fii.close()

#==================================================================
# Get SR84 data
#==================================================================

    if anl.addSR84E5:

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       files["SR84"],timestrs[c])

      # Get units
      SR84E5units = fh[0].variables["MSHR"].units

      # Get an area
      if anl.addfull:
        SR84E5 = E5fns.get_E5_ss_2D_var(
               fh,"MSHR",timi,loni,lati,times,ctime)
        SR84E5_nr = np.where(TCRWE5>0.001,np.nan,SR84E5)
        varnan = SR84E5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          SR84_mean_nr_E5[c] = np.nanmean(SR84E5_nr)
      if anl.addanom:
        SR84anomE5 = E5fns.get_E5_ss_2D_var(
               fh,"MSHR_anom",timi,loni,lati,times,ctime)
        SR84anomE5_nr = np.where(TCRWE5>0.001,np.nan,SR84anomE5)
        varnan = SR84anomE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          SR84anom_mean_nr_E5[c] = np.nanmean(SR84anomE5_nr)
      if anl.addanomn:
        SR84anomnE5 = E5fns.get_E5_ss_2D_var(
               fh,"MSHR_anomn",timi,loni,lati,times,ctime)
        SR84anomnE5_nr = np.where(TCRWE5>0.001,np.nan,SR84anomnE5)
        varnan = SR84anomnE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          SR84anomn_mean_nr_E5[c] = np.nanmean(SR84anomnE5_nr)

      # Close file
      for fii in fh:
        fii.close()

#==================================================================
# Get SR14 data
#==================================================================

    if anl.addSR14E5:

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       files["SR14"],timestrs[c])

      # Get units
      SR14E5units = fh[0].variables["MSHR"].units

      # Get an area
      if anl.addfull:
        SR14E5 = E5fns.get_E5_ss_2D_var(
               fh,"MSHR",timi,loni,lati,times,ctime)
        SR14E5_nr = np.where(TCRWE5>0.001,np.nan,SR14E5)
        varnan = SR14E5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          SR14_mean_nr_E5[c] = np.nanmean(SR14E5_nr)
      if anl.addanom:
        SR14anomE5 = E5fns.get_E5_ss_2D_var(
               fh,"MSHR_anom",timi,loni,lati,times,ctime)
        SR14anomE5_nr = np.where(TCRWE5>0.001,np.nan,SR14anomE5)
        varnan = SR14anomE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          SR14anom_mean_nr_E5[c] = np.nanmean(SR14anomE5_nr)
      if anl.addanomn:
        SR14anomnE5 = E5fns.get_E5_ss_2D_var(
               fh,"MSHR_anomn",timi,loni,lati,times,ctime)
        SR14anomnE5_nr = np.where(TCRWE5>0.001,np.nan,SR14anomnE5)
        varnan = SR14anomnE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          SR14anomn_mean_nr_E5[c] = np.nanmean(SR14anomnE5_nr)

      # Close file
      for fii in fh:
        fii.close()

#==================================================================
# Get CV18 data
#==================================================================

    if anl.addCV18E5:

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       files["CV18"],timestrs[c])

      # Get units
      CV18E5units = fh[0].variables["CONV"].units

      # Get an area
      if anl.addfull:
        CV18E5 = E5fns.get_E5_ss_2D_var(
               fh,"CONV",timi,loni,lati,times,ctime)
        CV18E5_nr = np.where(TCRWE5>0.001,np.nan,CV18E5)
        varnan = CV18E5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          CV18_mean_nr_E5[c] = np.nanmean(CV18E5_nr)
      if anl.addanom:
        CV18anomE5 = E5fns.get_E5_ss_2D_var(
               fh,"CONV_anom",timi,loni,lati,times,ctime)
        CV18anomE5_nr = np.where(TCRWE5>0.001,np.nan,CV18anomE5)
        varnan = CV18anomE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          CV18anom_mean_nr_E5[c] = np.nanmean(CV18anomE5_nr)
      if anl.addanomn:
        CV18anomnE5 = E5fns.get_E5_ss_2D_var(
               fh,"CONV_anomn",timi,loni,lati,times,ctime)
        CV18anomnE5_nr = np.where(TCRWE5>0.001,np.nan,CV18anomnE5)
        varnan = CV18anomnE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          CV18anomn_mean_nr_E5[c] = np.nanmean(CV18anomnE5_nr)

      # Close file
      for fii in fh:
        fii.close()

#==================================================================
# Get CV31 data
#==================================================================

    if anl.addCV31E5:

      # Find file and time indices
      fh,timi,times,ctime = E5fns.get_E5_ss_2D_fiti(
                       files["CV31"],timestrs[c])

      # Get units
      CV31E5units = fh[0].variables["CONV"].units

      # Get an area
      if anl.addfull:
        CV31E5 = E5fns.get_E5_ss_2D_var(
               fh,"CONV",timi,loni,lati,times,ctime)
        CV31E5_nr = np.where(TCRWE5>0.001,np.nan,CV31E5)
        varnan = CV31E5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          CV31_mean_nr_E5[c] = np.nanmean(CV31E5_nr)
      if anl.addanom:
        CV31anomE5 = E5fns.get_E5_ss_2D_var(
               fh,"CONV_anom",timi,loni,lati,times,ctime)
        CV31anomE5_nr = np.where(TCRWE5>0.001,np.nan,CV31anomE5)
        varnan = CV31anomE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          CV31anom_mean_nr_E5[c] = np.nanmean(CV31anomE5_nr)
      if anl.addanomn:
        CV31anomnE5 = E5fns.get_E5_ss_2D_var(
               fh,"CONV_anomn",timi,loni,lati,times,ctime)
        CV31anomnE5_nr = np.where(TCRWE5>0.001,np.nan,CV31anomnE5)
        varnan = CV31anomnE5_nr.flatten()
        if np.isnan(varnan).sum()/len(varnan)<anl.avgmissfrac:
          CV31anomn_mean_nr_E5[c] = np.nanmean(CV31anomnE5_nr)

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
    np.float64,timeunits,fileout,tE5)

#==================================================================
# Write TCWV information to file
#==================================================================

  if anl.addTCWVE5:

    description = "Total column water vapor from the ERA5 dataset averaged over a 5x5 degree area, centered on the precipitation system centroid, and areas with substantial rain water removed. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    if anl.addfull:
      mfns.write_var("TCWV_mean_nr_E5",
       "ERA5 Mean Total Column Water Vapor (Rain cells removed)",description,("tE5"),
       np.float64,TCWVE5units,fileout,TCWV_mean_nr_E5)
    if anl.addanom:
      mfns.write_var("TCWVanom_mean_nr_E5",
       "ERA5 Mean Total Column Water Vapor Anomaly (Rain cells removed)",description,("tE5"),
       np.float64,TCWVE5units,fileout,TCWVanom_mean_nr_E5)
    if anl.addanomn:
      mfns.write_var("TCWVanomn_mean_nr_E5",
       "ERA5 Mean Total Column Water Vapor Normalized Anomaly (Rain cells removed)",description,
       ("tE5"),np.float64,"",fileout,TCWVanomn_mean_nr_E5)

#==================================================================
# Write SH18 information to file
#==================================================================

  if anl.addSH18E5:

    description = "Specific humidity (1000-850hPa) from the ERA5 dataset averaged over a 5x5 degree area, centered on the precipitation system centroid, and areas with substantial rain water removed. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    if anl.addfull:
      mfns.write_var("SH18_mean_nr_E5",
       "ERA5 Mean Specific humidity (1000-850hPa) (Rain cells removed)",description,("tE5"),
       np.float64,SH18E5units,fileout,SH18_mean_nr_E5)
    if anl.addanom:
      mfns.write_var("SH18anom_mean_nr_E5",
       "ERA5 Mean Specific humidity (1000-850hPa) Anomaly (Rain cells removed)",description,("tE5"),
       np.float64,SH18E5units,fileout,SH18anom_mean_nr_E5)
    if anl.addanomn:
      mfns.write_var("SH18anomn_mean_nr_E5",
       "ERA5 Mean Specific humidity (1000-850hPa) Normalized Anomaly (Rain cells removed)",description,
       ("tE5"),np.float64,"",fileout,SH18anomn_mean_nr_E5)

#==================================================================
# Write SH84 information to file
#==================================================================

  if anl.addSH84E5:

    description = "Specific humidity (800-400hPa) from the ERA5 dataset averaged over a 5x5 degree area, centered on the precipitation system centroid, and areas with substantial rain water removed. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    if anl.addfull:
      mfns.write_var("SH84_mean_nr_E5",
       "ERA5 Mean Specific humidity (800-400hPa) (Rain cells removed)",description,("tE5"),
       np.float64,SH84E5units,fileout,SH84_mean_nr_E5)
    if anl.addanom:
      mfns.write_var("SH84anom_mean_nr_E5",
       "ERA5 Mean Specific humidity (800-400hPa) Anomaly (Rain cells removed)",description,("tE5"),
       np.float64,SH84E5units,fileout,SH84anom_mean_nr_E5)
    if anl.addanomn:
      mfns.write_var("SH84anomn_mean_nr_E5",
       "ERA5 Mean Specific humidity (800-400hPa) Normalized Anomaly (Rain cells removed)",description,
       ("tE5"),np.float64,"",fileout,SH84anomn_mean_nr_E5)

#==================================================================
# Write SR18 information to file
#==================================================================

  if anl.addSR18E5:

    description = "Shear magnitude (1000-850hPa) from the ERA5 dataset averaged over a 5x5 degree area, centered on the precipitation system centroid, and areas with substantial rain water removed. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    if anl.addfull:
      mfns.write_var("SR18_mean_nr_E5",
       "ERA5 Mean shear magnitude (1000-850hPa) (Rain cells removed)",description,("tE5"),
       np.float64,SR18E5units,fileout,SR18_mean_nr_E5)
    if anl.addanom:
      mfns.write_var("SR18anom_mean_nr_E5",
       "ERA5 Mean shear magnitude (1000-850hPa) Anomaly (Rain cells removed)",description,("tE5"),
       np.float64,SR18E5units,fileout,SR18anom_mean_nr_E5)
    if anl.addanomn:
      mfns.write_var("SR18anomn_mean_nr_E5",
       "ERA5 Mean shear magnitude (1000-850hPa) Normalized Anomaly (Rain cells removed)",description,
       ("tE5"),np.float64,"",fileout,SR18anomn_mean_nr_E5)

#==================================================================
# Write SR84 information to file
#==================================================================

  if anl.addSR84E5:

    description = "shear magnitude (800-400hPa) from the ERA5 dataset averaged over a 5x5 degree area, centered on the precipitation system centroid, and areas with substantial rain water removed. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    if anl.addfull:
      mfns.write_var("SR84_mean_nr_E5",
       "ERA5 Mean shear magnitude (800-400hPa) (Rain cells removed)",description,("tE5"),
       np.float64,SR84E5units,fileout,SR84_mean_nr_E5)
    if anl.addanom:
      mfns.write_var("SR84anom_mean_nr_E5",
       "ERA5 Mean shear magnitude (800-400hPa) Anomaly (Rain cells removed)",description,("tE5"),
       np.float64,SR84E5units,fileout,SR84anom_mean_nr_E5)
    if anl.addanomn:
      mfns.write_var("SR84anomn_mean_nr_E5",
       "ERA5 Mean shear magnitude (800-400hPa) Normalized Anomaly (Rain cells removed)",description,
       ("tE5"),np.float64,"",fileout,SR84anomn_mean_nr_E5)

#==================================================================
# Write SR14 information to file
#==================================================================

  if anl.addSR14E5:

    description = "shear magnitude (800-400hPa) from the ERA5 dataset averaged over a 5x5 degree area, centered on the precipitation system centroid, and areas with substantial rain water removed. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    if anl.addfull:
      mfns.write_var("SR14_mean_nr_E5",
       "ERA5 Mean shear magnitude (800-400hPa) (Rain cells removed)",description,("tE5"),
       np.float64,SR14E5units,fileout,SR14_mean_nr_E5)
    if anl.addanom:
      mfns.write_var("SR14anom_mean_nr_E5",
       "ERA5 Mean shear magnitude (800-400hPa) Anomaly (Rain cells removed)",description,("tE5"),
       np.float64,SR14E5units,fileout,SR14anom_mean_nr_E5)
    if anl.addanomn:
      mfns.write_var("SR14anomn_mean_nr_E5",
       "ERA5 Mean shear magnitude (800-400hPa) Normalized Anomaly (Rain cells removed)",description,
       ("tE5"),np.float64,"",fileout,SR14anomn_mean_nr_E5)

#==================================================================
# Write CV18 information to file
#==================================================================

  if anl.addCV18E5:

    description = "convergence (1000-850hPa) from the ERA5 dataset averaged over a 5x5 degree area, centered on the precipitation system centroid, and areas with substantial rain water removed. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    if anl.addfull:
      mfns.write_var("CONV18_mean_nr_E5",
       "ERA5 Mean convergence (1000-850hPa) (Rain cells removed)",description,("tE5"),
       np.float64,CV18E5units,fileout,CV18_mean_nr_E5)
    if anl.addanom:
      mfns.write_var("CONV18anom_mean_nr_E5",
       "ERA5 Mean convergence (1000-850hPa) Anomaly (Rain cells removed)",description,("tE5"),
       np.float64,CV18E5units,fileout,CV18anom_mean_nr_E5)
    if anl.addanomn:
      mfns.write_var("CONV18anomn_mean_nr_E5",
       "ERA5 Mean convergence (1000-850hPa) Normalized Anomaly (Rain cells removed)",description,
       ("tE5"),np.float64,"",fileout,CV18anomn_mean_nr_E5)

#==================================================================
# Write CV31 information to file
#==================================================================

  if anl.addCV31E5:

    description = "Convergence (300-100hPa) from the ERA5 dataset averaged over a 5x5 degree area, centered on the precipitation system centroid, and areas with substantial rain water removed. Exact latitude and longitude coordinates are in the attributes of the groups lonsE5 and latsE5."
    if anl.addfull:
      mfns.write_var("CONV18_mean_nr_E5",
       "ERA5 Mean Convergence (300-100hPa) (Rain cells removed)",description,("tE5"),
       np.float64,CV31E5units,fileout,CV31_mean_nr_E5)
    if anl.addanom:
      mfns.write_var("CONV31anom_mean_nr_E5",
       "ERA5 Mean Convergence (300-100hPa) Anomaly (Rain cells removed)",description,("tE5"),
       np.float64,CV31E5units,fileout,CV31anom_mean_nr_E5)
    if anl.addanomn:
      mfns.write_var("CONV31anomn_mean_nr_E5",
       "ERA5 Mean Convergence (300-100hPa) Normalized Anomaly (Rain cells removed)",description,
       ("tE5"),np.float64,"",fileout,CV31anomn_mean_nr_E5)

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
  print("Begin serial loop over TIPS")
  for fn in range(len(filenamesrun)):
    driver_addvars(fn)

# Parrallel loop over PFs
if anl.serialorparallel==2:
  print("Begin parrallel loop over TIPS")
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
