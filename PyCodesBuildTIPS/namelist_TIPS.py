#==================================================================
# Namelist for running TIPS programs
#
# In here you will input all namelist options for the 
#  python programs involved in :
#  * create_FiT_input_files
#  * process_FiTobs
#  * add_vars
#  * add_ERA5_vars
#
# Namelist must be in same directory as python scripts
#
# James Russell 2020
#==================================================================

# Import Python libraries (don't change)
import numpy as np

#==================================================================
# General Namelist
#==================================================================
# Do not change
def general():
  pass
#==================================================================

# Directory for custom functions
general.fnsdir = "/uufs/chpc.utah.edu/common/home/u0816744/general_functions/"

# Directory and filename for input IMERG data
general.datadirin = "/uufs/chpc.utah.edu/common/home/varble-group2/IMERG/"
general.fileidin  = "3B-HHR.MS.MRG.3IMERG."
general.datahdf5  = True
general.datanc4   = False

# Directory and filename for TIPS files
general.datadirTIPS = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_CPEX-AW/TIPS_test/"
general.fileidTIPS = "TIPS_"

#==================================================================
# Namelist for create_FiT_input_files
#==================================================================
# Do not change
def create():
  pass
#==================================================================

# Date and time range
create.starttime = "20180601"
create.endtime   = "20180610" # Actual last day is day before

# Thresholds
create.tholdtype = 2 # Threshold type (1 = fixed, 2 = normalized)
create.tholds = [0.25,0.5] # Thresholds
create.minthold = 1. # Min precip threshold. Only required for tholdtype=2.

# Directory and filename for output FiT input data
create.datadirout = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_CPEX-AW/FiT_input_test/"
create.fileidout  = "IMERG_FiT_tholds_"

# Subset regions (ranges have no affect if ssreg=False)
create.ssreg  = True
create.latN   = 30
create.latS   = -10
create.lonW   = -70
create.lonE   = 0
create.ssname = "Atl"

# Smoothing
create.smooth   = True
create.gaussian = False
create.stddev   = 1.5  # Standard deviation for gaussian smoothing
create.uniform  = True
create.width    = 5    # Number of points to spread running average

# Parallelization
create.serialorparallel = 2 # serial=1 parallel=2
create.njobs = 32 # Number of cores for parallelization

#==================================================================
# Namelist for process_FiTobs
#==================================================================
# Do not change
def process():
  pass
#==================================================================

# Directory and filename for FiT input files
# Should obtain from create_FiT_input_files.py
process.datadirFin = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_CPEX-AW/FiT_input_test/"
process.fileidFin  = "IMERG_FiT_tholds_Atl_"

# Directory and filename for FiT output data files
process.datadirFi  = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_CPEX-AW/FiT_output_test/"
process.fileidFi1  = "IMERG_tracked_"
process.fileidFi2  = "_4Dobjects.nc"
process.fileidtxt  = "4Dobject_tree.txt"

# Only process a set of objects
process.subsetobs  = True
process.ob1        = 100000 # First obid
process.ob2        = 100100 # Last obid
process.subsetsztm = True
process.nsz        = 30     # Minimum no of pixels (30 advised)
process.nt         = 6      # Minimum number of times

# Specify merging distance from tracking (number of pixels)
process.mergedist = 30

# Specify reference time in format YYYY-MM-DD hh:mm:ss
process.reftime = "1900-01-01 00:00:00"

# Parallelization
process.serialorparallel = 2 # 1=serial; 2=parallel
process.njobs = 32 # Number of jobs for parallel

#==================================================================
# Namelist for add_vars
#==================================================================
# Do not change
def addvars():
  pass
#==================================================================

# Subset (for certain range of dates) 
addvars.ssdat = False
addvars.date1 = "20180601"
addvars.date2 = "20180602"
addvars.ssobs = False
addvars.obid1 = "100000"
addvars.obid2 = "100010"

# Number of processes for parrallelization
addvars.serialorparallel = 2
addvars.njobs = 32

# Variables desired
addvars.addmaxrr          = True # Maximum rain rate
addvars.addmeanrr         = True # Mean rain rate
addvars.addmedianrr       = True # Median rain rate
addvars.addstddevrr       = True # Standard deviation of the rain
                           #  rates
addvars.addpieces         = True # Number of pieces 
addvars.addpiecesc        = True # As above but for convective rain pixels
addvars.addarea           = True # Area of the PF
addvars.addvrr            = True # Volumetric rain rate
addvars.addpropagation    = True # Propagation characteristics
addvars.addTCinfo         = True # Flags indicating proximity to 
                           #  TC center
addvars.addlandinfo       = True # Flags indicating locations over 
                           #  land
addvars.addboundaryinfo   = True # Time-series indicating if PF 
                           #  touching domain boundary
addvars.addlocaltime      = True # Local solar time of the PF
addvars.addasymmetry      = True # Asymmetry shape parameter 
                           #  (Zick et al. 2016)
addvars.addasymmetryc     = True # As above for convective pixels
addvars.addfragmentation  = True # Fragmentation shape parameter 
                           #  (Zick et al. 2016)
addvars.addfragmentationc = True # As above for convective pixels
addvars.adddispersion     = True # Dispersion shape parameter 
                           #  (Zick et al. 2016)
addvars.adddispersionc    = True # As above for convective pixels
addvars.addaxesshape      = True # Array of variables based on 
                           #  the major and minor axes from 
                           #  eigenvalue/vectors
addvars.addaxesshapec     = True # As above for convective pixels
addvars.addperimeter      = True # Distance around perimeter of 
                           #  shape (alpha-shape method)
addvars.addconvrain       = True # Flags indicating whether rain 
                           #  is convective
addvars.addconvarea       = True # Area of the convective region 
                           #  (addconvrain must also be True)
addvars.addconvvrr        = True # Volume of convective rainfall 
                           #  (addconvrain must also be True)

# Inputs for specific variables
  
# Grid spacing in degrees lon, lat (only required for 
#  area/volrainrate/shape)
addvars.dx = 0.1
addvars.dy = 0.1

# Directory and filename of TC data (only required for TC 
#  information)
addvars.dataTCdir = "/uufs/chpc.utah.edu/common/home/varble-group2/IBTrACS/"
addvars.fileTCid  = "IBTrACS.ALL.v04r00.nc"

# Convective rain rate threshold (mm/hr)
addvars.convrainthold = 10
  
# Shape metric indexes
addvars.minshapesize  = 30  # Minimum number of pixels for PF shape
addvars.minshapefrag  = 0.3 # Minimum fragmentation for PF shape
addvars.minshapedisp  = 0.5 # Minimum dispersion for PF shape
addvars.minshapesizec = 30  # Minimum number of pixels for PF shape 
addvars.minshapefragc = 0.5 # Minimum fragmentation for convective shape
addvars.minshapedispc = 0.8 # Minimum dispersion for convective shape

#==================================================================
# Namelist for add_ERA5_vars
#==================================================================
# Do not change
def addERA5vars():
  pass
#==================================================================

# Subset (for certain range of dates) 
addERA5vars.ssdat = False
addERA5vars.date1 = "20180601"
addERA5vars.date2 = "20180602"
addERA5vars.ssobs = False
addERA5vars.obid1 = "103231"
addERA5vars.obid2 = "103233"

# Number of processes for parrallelization
addERA5vars.serialorparallel = 2
addERA5vars.njobs = 20

# Type of area or mean
addERA5vars.addctarea = True # Area centered on TIPS
addERA5vars.addctmean = True  # Mean of centered area
addERA5vars.addinarea = False # All points in calculated inflow
addERA5vars.addinmean = False # Add mean of inflow region

# Rain check or no rain check
addERA5vars.addrainchk   = True
addERA5vars.addnorainchk = True

# Variables desired
addERA5vars.addTCWVE5    = True # ERA5 Total Column Water Vapor
addERA5vars.addCAPEE5    = True # ERA5 Convective Available Potential Energy
addERA5vars.addSR18E5    = True # ERA5 Boundary Layer Shear
addERA5vars.addSR82E5    = True # ERA5 Free-troposphere Shear
addERA5vars.addSH18E5    = True # ERA5 Boundary Layer Specific Humidity
addERA5vars.addSH82E5    = True # ERA5 Free-troposphere Specific Humidity

# ERA5 domain variables
addERA5vars.hda           = 2.5 # Half data area in degrees
addERA5vars.hoursbefore   = np.arange(48,0,-3) # Hours before (descending)
addERA5vars.hoursafter    = addERA5vars.hoursbefore[::-1] # Hours after (ascending)
addERA5vars.avgmissfrac   = 0.8 # fraction of domain missing before average is not carried out

# Directory and filenames of ERA5 data
addERA5vars.dataE5dir     = "/uufs/chpc.utah.edu/common/home/varble-group1/ERA5/"
addERA5vars.fileTCWVE5id  = "moisture/ERA5.TCWV."
addERA5vars.fileCAPEE5id  = "convparams/ERA5.CAPE."
addERA5vars.fileUSR18E5id = "shear/ERA5.USHR_1000-850hPamean."
addERA5vars.fileVSR18E5id = "shear/ERA5.VSHR_1000-850hPamean."
addERA5vars.fileUSR82E5id = "shear/ERA5.USHR_850-200hPamean."
addERA5vars.fileVSR82E5id = "shear/ERA5.VSHR_850-200hPamean."
addERA5vars.fileSH18E5id  = "moisture/ERA5.SPHU_1000-850hPamean."
addERA5vars.fileSH82E5id  = "moisture/ERA5.SPHU_850-200hPamean."
addERA5vars.fileTCRWE5id  = "clouds/ERA5.TCRW."

#==================================================================
# End namelist
#==================================================================

