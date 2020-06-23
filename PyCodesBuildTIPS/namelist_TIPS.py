#==========================================================
# Namelist for running TIPS programs
# 
# Must be in directory 
#
# In here you will input all namelist options for the 
#  python programs involved in the FiT:
#  * create_FiT_input_files
#  * process_FiTobs
#  * add_vars
#  * add_ERA5_vars
#==========================================================

# Import Python libraries (don't change)
import numpy as np

#==========================================================
# General Namelist
#==========================================================

# Directory for custom functions
fnsdir = "/uufs/chpc.utah.edu/common/home/u0816744/general_functions/"

# Directory and filename for input IMERG data
datadirin = "/uufs/chpc.utah.edu/common/home/varble-group2/IMERG/"
fileidin  = "3B-HHR.MS.MRG.3IMERG."
datahdf5  = True
datanc4   = False

# Directory and filename for TIPS files
datadirTIPS = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_CPEX-AW/TIPS_test/"
fileidTIPS = "TIPS_"

#==========================================================
# Namelist for create_FiT_input_files
#==========================================================

# Date and time range
starttime = "20180601"
endtime   = "20180610" # Actual last day is day before

# Thresholds
tholdtype = 2 # Threshold type (1 = fixed, 2 = normalized)
tholds = [0.25,0.5] # Thresholds
minthold = 1. # Min precip threshold. Only required for tholdtype=2.

# Directory and filename for output FiT input data
datadiroutc = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_CPEX-AW/FiT_input_test/"
fileidoutc  = "IMERG_FiT_tholds_"

# Subset regions (ranges have no affect if ssreg=False)
ssreg  = True
latN   = 30
latS   = -10
lonW   = -70
lonE   = 0
ssname = "Atl"

# Smoothing
smooth   = True
gaussian = False
stddev   = 1.5  # Standard deviation for gaussian smoothing
uniform  = True
width    = 5    # Number of points to spread running average

# Parallelization
serialorparallelc = 2 # serial=1 parallel=2
njobsc = 32 # Number of cores for parallelization

#============================================================
# Namelist for process_FiTobs
#============================================================

# Directory and filename for FiT input files
# Should obtain from create_FiT_input_files.py
datadirFin = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_CPEX-AW/FiT_input_test/"
fileidFin  = "IMERG_FiT_tholds_Atl_"

# Directory and filename for FiT output data files
datadirFi  = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_CPEX-AW/FiT_output_test/"
fileidFi1  = "IMERG_tracked_"
fileidFi2  = "_4Dobjects.nc"
fileidtxt  = "4Dobject_tree.txt"

# Only process a set of objects
subsetobs = True
ob1       = 100000 # First obid
ob2       = 100100 # Last obid
subsetsztm= True
nsz       = 30     # Minimum no of pixels (30 advised)
nt        = 6      # Minimum number of times

# Specify merging distance from tracking (number of pixels)
mergedist = 30

# Specify reference time in format YYYY-MM-DD hh:mm:ss
reftime = "1900-01-01 00:00:00"

# Number of processes for parallelization
serialorparallelp = 2
njobsp = 32

#==================================================================
# Namelist for add_vars
#==================================================================

# Subset (for certain range of dates) 
ssdat = False
date1 = "20180601"
date2 = "20180602"
ssobs = False
obid1 = "100000"
obid2 = "100010"

# Number of processes for parrallelization
serialorparallela = 2
njobsa = 32

# Variables desired
addmaxrr          = True # Maximum rain rate
addmeanrr         = True # Mean rain rate
addmedianrr       = True # Median rain rate
addstddevrr       = True # Standard deviation of the rain
                         #  rates
addpieces         = True # Number of pieces 
addpiecesc        = True # As above but for convective rain pixels
addarea           = True # Area of the PF
addvrr            = True # Volumetric rain rate
addpropagation    = True # Propagation characteristics
addTCinfo         = True # Flags indicating proximity to 
                         #  TC center
addlandinfo       = True # Flags indicating locations over 
                         #  land
addboundaryinfo   = True # Time-series indicating if PF 
                         #  touching domain boundary
addlocaltime      = True # Local solar time of the PF
addasymmetry      = True # Asymmetry shape parameter 
                         #  (Zick et al. 2016)
addasymmetryc     = True # As above for convective pixels
addfragmentation  = True # Fragmentation shape parameter 
                         #  (Zick et al. 2016)
addfragmentationc = True # As above for convective pixels
adddispersion     = True # Dispersion shape parameter 
                         #  (Zick et al. 2016)
adddispersionc    = True # As above for convective pixels
addaxesshape      = True # Array of variables based on 
                         #  the major and minor axes from 
                         #  eigenvalue/vectors
addaxesshapec     = True # As above for convective pixels
addperimeter      = True # Distance around perimeter of 
                         #  shape (alpha-shape method)
addconvrain       = True # Flags indicating whether rain 
                         #  is convective
addconvarea       = True # Area of the convective region 
                         #  (addconvrain must also be True)
addconvvrr        = True # Volume of convective rainfall 
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
minshapesize  = 30  # Minimum number of pixels for PF shape
minshapefrag  = 0.3 # Minimum fragmentation for PF shape
minshapedisp  = 0.5 # Minimum dispersion for PF shape
minshapesizec = 30  # Minimum number of pixels for PF shape
minshapefragc = 0.5 # Minimum fragmentation for convective shape
minshapedispc = 0.8 # Minimum dispersion for convective shape

#==================================================================
# Namelist for add_ERA5_vars
#==================================================================

# Subset (for certain range of dates) 
ssdataE = False
date1aE = "20180601"
date2aE = "20180602"
ssobsaE = False
obid1aE = "103231"
obid2aE = "103233"

# Number of processes for parrallelization
serialorparallelaE = 2
njobsaE = 20

# Type of area or mean
addctarea = True # Area centered on TIPS
addctmean = True  # Mean of centered area
addinarea = False # All points in calculated inflow
addinmean = False # Add mean of inflow region

# Rain check or no rain check
addrainchk   = True
addnorainchk = True

# Variables desired
addTCWVE5    = True # ERA5 Total Column Water Vapor
addCAPEE5    = True # ERA5 Convective Available Potential Energy
addSR18E5    = True # ERA5 Boundary Layer Shear
addSR82E5    = True # ERA5 Free-troposphere Shear
addSH18E5    = True # ERA5 Boundary Layer Specific Humidity
addSH82E5    = True # ERA5 Free-troposphere Specific Humidity

# ERA5 domain variables
hda           = 2.5 # Half data area in degrees
hoursbefore   = np.arange(48,0,-3) # Hours before (descending)
hoursafter    = hoursbefore[::-1] # Hours after (ascending)
avgmissfrac   = 0.8 # fraction of domain missing before average is not carried out

# Directory and filenames of ERA5 data
dataE5dir     = "/uufs/chpc.utah.edu/common/home/varble-group1/ERA5/"
fileTCWVE5id  = "moisture/ERA5.TCWV."
fileCAPEE5id  = "convparams/ERA5.CAPE."
fileUSR18E5id = "shear/ERA5.USHR_1000-850hPamean."
fileVSR18E5id = "shear/ERA5.VSHR_1000-850hPamean."
fileUSR82E5id = "shear/ERA5.USHR_850-200hPamean."
fileVSR82E5id = "shear/ERA5.VSHR_850-200hPamean."
fileSH18E5id  = "moisture/ERA5.SPHU_1000-850hPamean."
fileSH82E5id  = "moisture/ERA5.SPHU_850-200hPamean."
fileTCRWE5id  = "clouds/ERA5.TCRW."
