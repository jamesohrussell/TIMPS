#==================================================================
# Add variables to IMERG PF netcdf files.
#
# Must provide in namelist:
#  * Directory and filename identifiers for PF files
#  * Directory and filename identifiers for IBtracs data
#     (only if TC information is desired)
#  * Variables desired
#  * Grid spacing of input data pixels in degrees 
#     (IMERG = 0.1, only required if area or volumetric 
#     rain rate is desired)
#  * Subset information if desired
#
# James Russell 2019
#==================================================================

# Import python libraries (do not change)
from netCDF4 import Dataset
import glob
import datetime as dt
import time as tm
import driver_addvars as da
from joblib import Parallel, delayed
import os

#==================================================================
# Namelist
#==================================================================

# Directory for custom functions
fnsdir = "/uufs/chpc.utah.edu/common/home/u0816744/general_functions/"

# Directory and filename for TIPS files
datadir = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_CPEX-AW/TIPS_test/"
fileid  = "TIPS_"

# Subset (for certain range of dates) 
ssdat = False
date1 = "20180601"
date2 = "20180602"
ssobs = False
obid1 = "100000"
obid2 = "100010"

# Number of processes for parrallelization
serialorparallel = 1
njobs = 16

# Variables desired
addmaxrr          = False # Maximum rain rate
addmeanrr         = False # Mean rain rate
addmedianrr       = False # Median rain rate
addstddevrr       = False # Standard deviation of the rain
                          #  rates
addpieces         = True # Number of pieces 
addpiecesc        = False # As above but for convective rain pixels
addarea           = False # Area of the PF
addvrr            = False # Volumetric rain rate
addpropagation    = False # Propagation characteristics
addTCinfo         = False # Flags indicating proximity to 
                          #  TC center
addlandinfo       = False # Flags indicating locations over 
                          #  land
addboundaryinfo   = False # Time-series indicating if PF 
                          #  touching domain boundary
addlocaltime      = False # Local solar time of the PF
addasymmetry      = False # Asymmetry shape parameter 
                          #  (Zick et al. 2016)
addasymmetryc     = False # As above for convective pixels
addfragmentation  = False # Fragmentation shape parameter 
                          #  (Zick et al. 2016)
addfragmentationc = False # As above for convective pixels
adddispersion     = False # Dispersion shape parameter 
                          #  (Zick et al. 2016)
adddispersionc    = True # As above for convective pixels
addaxesshape      = True # Array of variables based on 
                          #  the major and minor axes from 
                          #  eigenvalue/vectors
addaxesshapec     = True # As above for convective pixels
addperimeter      = False # Distance around perimeter of 
                          #  shape (alpha-shape method)
addconvrain       = False # Flags indicating whether rain 
                          #  is convective
addconvarea       = False # Area of the convective region 
                          #  (addconvrain must also be True)
addconvvrr        = False # Volume of convective rainfall 
                          #  (addconvrain must also be True)

# Inputs for specific variables

# Grid spacing in degrees lon, lat (only required for 
#  area/volrainrate/shape)
dx = 0.1
dy = 0.1

# Directory and filename of TC data (only required for TC 
#  information)
dataTCdir = "/uufs/chpc.utah.edu/common/home/varble-group2/IBTrACS/"
fileTCid  = "IBTrACS.ALL.v04r00.nc"

# Convective rain rate threshold (mm/hr)
convrainthold = 10

# Shape metric indexes
minshapesize  = 30   # Minimum number of pixels for PF shape
minshapefrag  = 0.35 # Minimum fragmentation for PF shape
minshapedisp  = 0.7  # Minimum dispersion for PF shape
minshapesizec = 30   # Minimum number of pixels for PF shape
minshapefragc = 0.5 # Minimum fragmentation for convective shape
minshapedispc = 1.0  # Minimum dispersion for convective shape

#==================================================================
# Initialize timer
#==================================================================

startnow = tm.time()

#==================================================================
# Write namelist to a dictionary
#==================================================================

# Put namelist information in dictionaries
namelist = {}
namelist["fnsdir"] = str(fnsdir)

# Add which variables are selected
namelist["addmaxrr"] = str(addmaxrr)
namelist["addmeanrr"] = str(addmeanrr)
namelist["addmedianrr"] = str(addmedianrr)
namelist["addstddevrr"] = str(addstddevrr)
namelist["addpieces"] = str(addpieces)
namelist["addpiecesc"] = str(addpiecesc)
namelist["addarea"] = str(addarea)
namelist["addvrr"] = str(addvrr)
namelist["addpropagation"] = str(addpropagation)
namelist["addTCinfo"] = str(addTCinfo)
namelist["addlandinfo"] = str(addlandinfo)
namelist["addboundaryinfo"] = str(addboundaryinfo)
namelist["addlocaltime"] = str(addlocaltime)
namelist["addconvrain"] = str(addconvrain)
namelist["addconvarea"] = str(addconvarea)
namelist["addconvvrr"] = str(addconvvrr)
namelist["addaxesshape"] = str(addaxesshape)
namelist["addaxesshapec"] = str(addaxesshapec)
namelist["addperimeter"] = str(addperimeter)
namelist["addasymmetry"] = str(addasymmetry)
namelist["addasymmetryc"] = str(addasymmetryc)
namelist["addfragmentation"] = str(addfragmentation)
namelist["addfragmentationc"] = str(addfragmentationc)
namelist["adddispersion"] = str(adddispersion)
namelist["adddispersionc"] = str(adddispersionc)

if addarea or addvrr or addconvarea or addconvvrr or \
   addperimeter or addasymmetry or addfragmentation or \
   addaxesshape or addboundaryinfo or addasymmetryc or \
   addaxesshapec or addfragmentationc or adddispersion or \
   adddispersionc or addpieces:
  namelist["dx"] = dx
  namelist["dy"] = dy
if addasymmetry or addfragmentation or adddispersion or \
   addaxeshape:
  namelist["minshapesize"] = minshapesize
if addasymmetryc or addfragmentationc or adddispersionc or \
   addaxeshapec:
  namelist["minshapesizec"] = minshapesizec
if addfragmentation or adddispersion:
  namelist["minshapefrag"] = minshapefrag
  namelist["minshapedisp"] = minshapedisp
if addfragmentationc or adddispersionc:
  namelist["minshapefragc"] = minshapefragc
  namelist["minshapedispc"] = minshapedispc
if addTCinfo:
  namelist["dataTCdir"] = str(dataTCdir)
  namelist["fileTCid"] = str(fileTCid)
if addconvrain or addconvarea or addconvvrr or addasymmetryc or \
   addfragmentationc or addaxesshapec or adddispersionc:
  namelist["convrainthold"] = convrainthold

# Write namelist dictionary to netcdf file for reading 
#  during parallel loop
print("Writing namelist to netcdf file")
nlfileout = Dataset("namelist_av.nc","w",format="NETCDF4")
for k,v in namelist.items():
  setattr(nlfileout, k,  v)
nlfileout.close()

#==================================================================
# Generate list of IPF files to process
#==================================================================

# Reads directory and file names
if ssdat:

  print("Subsetting by date")

  # Generate a list of filenames with dates to search for
  print("Generating filenames to search for")
  start = dt.datetime.strptime(date1,"%Y%m%d")
  end = dt.datetime.strptime(date2,"%Y%m%d")
  datearr = (start + dt.timedelta(days=x) for x in range(0,(end-start).days))
  filen = []
  for dateobj in datearr:
    filen.append(datadir+fileid+"*"+dateobj.strftime("%Y%m%d")+"*")

  # Reads directory and filenames
  print("Getting directory and filenames")
  filenamesrun = []
  for f in filen:
      filenamesrun = filenamesrun+sorted(glob.glob(f))

elif ssobs:

  print("Subsetting by object ids")

  print("Generating filenames to search for")
  obids = [i for i in range(int(obid1),int(obid2))]
  filen = []
  for i in range(len(obids)):
    filen.append(datadir+fileid+str(obids[i])+"*")

  # Reads directory and filenames
  print("Getting directory and filenames")
  filenamesrun = []
  for f in filen:
    print("Adding "+str(f))
    filenamesrun = filenamesrun+sorted(glob.glob(f))

  #import os
  #for f in filen:
  #  print(os.system("ls "+str(f)))

else:

  print("No subsetting")

  # Find all files in directory
  print("Generating list of files")
  filenamesrun = sorted(glob.glob(datadir+fileid+'*'))

# Write file names to a text file for reading during 
#  parallel loop
print("Writing filenames to text file")
ffn = open("filenames_av.txt","w")
for i in range(len(filenamesrun)):
  ffn.write(filenamesrun[i]+"\n")
ffn.close()

#==================================================================
# Loop over IPF files in parrallel
#==================================================================

if serialorparallel==1:
  print("Begin serial loop over objects")
  for fn in range(len(filenamesrun)):
    da.driver_addvars(fn)

# Parrallel loop over PFs
if serialorparallel==2:
  print("Begin parrallel loop over IPFs")
  Parallel(n_jobs=njobs)(delayed(da.driver_addvars)(fn) \
    for fn in range(len(filenamesrun)))

#==================================================================
# Final clean up tasks
#==================================================================

# Remove namelist and filename files
print("Cleaning up")
os.system("rm filenames_av.txt")
os.system("rm namelist_av.nc")
os.system("rm -rf __pycache__")

# End timer
endnow = tm.time()
print("Program took "+str(endnow - startnow)+"s")

#==================================================================
# End code
#==================================================================
