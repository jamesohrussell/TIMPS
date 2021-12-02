#==================================================================
# Namelist for running TIMPS programs
#
# In here you will input all namelist options for the 
#  python programs and the tracking:
#  * create_FiT_input_files
#  * FiT tracking
#  * write_tracking_output
#  * process_FiTobs
#  * add_vars
#  * add_ERA5_vars
#
# Namelist must be in same directory as python scripts
#
# James Russell 2021
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
general.datadirin = "/uufs/chpc.utah.edu/common/home/\
varble-group2/IMERG/final/"
general.fileidin  = "3B-HHR.MS.MRG.3IMERG."
general.datahdf5  = True
general.datanc4   = False

# Directory and filenames
general.datadirFiTin  = \
"/uufs/chpc.utah.edu/common/home/zipser-group2/TIMPS/processing_data_and_scripts/FiTin_2012/"
general.fileidFiTin   = "IMERG_FiT_tholds_"
general.datadirFiTout = \
"/uufs/chpc.utah.edu/common/home/zipser-group2/TIMPS/processing_data_and_scripts/FiTout_2012/"
general.fileidFiTout1 = "IMERG_tracked_"
general.fileidFiTout2 = "_4Dobjects.nc"
general.fileidFiToutT = "4Dobject_tree.txt"
general.datadirtrkout = \
"/uufs/chpc.utah.edu/common/home/zipser-group2/TIMPS/processing_data_and_scripts/Tracking_2012/"
general.fileidtrkout  = "Tracking_"
general.datadirTIMPS   = \
"/uufs/chpc.utah.edu/common/home/zipser-group2/TIMPS/TIMPS_files/2012/"
general.fileidTIMPS = "TIMPS_"

# Date and time range
general.starttime = "20111225" # Needs a spin-up period (i.e. 7 days)
general.endtime   = "20130115" # Actual last date is day before

# Split distance: number of pixels centers of two objects that 
#  were formerly one can be apart before splitting
general.splitdist = 20

# Is domain periodic in zonal direction
general.domperiodicx = True

# Convective rain rate threshold (mm/hr)
general.convrainthold = 10

# Parallelization for writing tracking output
general.wtserialorparallel = 2 # serial=1 parallel=2
general.wtnjobs = 4

#==================================================================
# Namelist for create_FiT_input_files
#==================================================================
# Do not change
def create():
  pass
#==================================================================

# Thresholds
# Threshold type (1 = multiple fixed; 2 = normalized; 
#  3 = contiguous area)
create.tholdtype = 2
# Thresholds. Required for tholdtype=1,2.
create.tholds = [0.11,0.33]
# Min precip threshold. Required for tholdtype=2,3.
create.minthold = 1.

# Subset regions (ranges have no affect if ssreg=False)
create.ssreg  = True
create.latN   = 30.05
create.latS   = -30.05
create.lonW   = -180
create.lonE   = 180
create.ssname = "GlobalTropics"

# Smoothing
create.smooth   = True
create.gaussian = False
create.stddev   = 1.5  # Standard deviation for gaussian smoothing
create.uniform  = True
create.width    = 5    # Number of points to spread running average

# Parallelization
# Parallization benchmarks (time for running 10 days of tracking)
# 52 = 10.75s
# 32 = 10.2s
# 20 = 10.53s
# 16 = 10.6s
# 10 = 11.98s
# 8  = 13.24s
create.serialorparallel = 2 # serial=1 parallel=2
create.njobs = 16 # Number of cores for parallelization

#==================================================================
# Namelist for process_FiTobs
#==================================================================
# Do not change
def process():
  pass
#==================================================================

# Only process a set of objects
process.subsetobs = False
process.ob1       = 240427 # First obid
process.ob2       = 240430 # Last obid
process.subsetsz  = True
process.nsz       = 30     # Minimum no of pixels (30 advised)
process.subsettm  = True
process.nt        = 12     # Minimum number of times
process.subsetmrn = True   # Minimum rain threshold (set in general with convective rain threshold)
process.subsetdts = True
process.date1     = "201201010000" #YYYYMMDDhhmm
process.date2     = "201212312330"

# Specify reference time in format YYYY-MM-DD hh:mm:ss
process.reftime = "1900-01-01 00:00:00"

# Parallelization
# Parallization benchmarks (time for running 2 weeks of global tropical strip tracking)
# 64 = 86s
# 20 = 27m54s
# 10 = 
# serial = 25m43s
process.serialorparallel = 1 # 1=serial; 2=parallel
process.njobs = 4 # Number of jobs for parallel

#==================================================================
# Namelist for add_vars
#==================================================================
# Do not change
def addvars():
  pass
#==================================================================

# Subset (for certain range of dates) 
addvars.ssdat = False
addvars.date1 = "20120101"
addvars.date2 = "20130131"
addvars.ssobs = False
addvars.obid1 = "1000000"
addvars.obid2 = "1000200"
addvars.sslst = False
addvars.lstfn = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/james/data/UG_file_lists/fn_ug_all_80_2012.txt"

# Number of processes for parrallelization
addvars.serialorparallel = 1
addvars.njobs = 8

# Variables desired
addvars.addlocaltime      = True# Local solar time of the PF
addvars.addmaxrr          = True# Maximum rain rate
addvars.addmeanrr         = True# Mean of non-zero rain rates
addvars.addmeanrr1mmhr    = True# Mean of >1mm/hr rain rates
addvars.addmeanrr10mmhr   = True# Mean of >10mm/hr rain rates
addvars.addmedianrr       = True# Median non-zero rain rate
addvars.addmedianrr1mmhr  = True# Median of >1mm/hr rain rates
addvars.addmedianrr10mmhr = True# Median of >10mm/hr rain rates
addvars.addstddevrr       = True# Standard deviation of non-zero rain rates
addvars.addstddevrr1mmhr  = True# Standard deviation of >1mm/hr rain rates
addvars.addstddevrr10mmhr = True# Standard deviation of >10mm/hr rain rates
addvars.addskewrr         = True# Skewness of non-zero rain rates
addvars.addskewrr1mmhr    = True# Skewness of >1mm/hr rain rates
addvars.addskewrr10mmhr   = True# Skewness of >10mm/hr rain rates
addvars.addpieces         = True# Number of pieces non-zero rain rates
addvars.addpieces1mmhr    = True# Number of pieces >1mm/hr rain rates
addvars.addpieces10mmhr   = True# Number of pieces >10mm/hr rain rates
addvars.addarea           = True# Area of non-zero rain rates
addvars.addarea1mmhr      = True# Area of >1mm/hr rain rates
addvars.addarea10mmhr     = True# Area of >10mm/hr rain rates
addvars.addvrr            = True# Volumetric of non-zero rain rates
addvars.addvrr1mmhr       = True# Volumetric of >1mm/hr rain rates
addvars.addvrr10mmhr      = True# Volumetric of >10mm/hr rain rates
addvars.addpropagation    = True# Propagation characteristics
addvars.addlandinfo       = True# Flags indicating locations over land
addvars.addTCinfo         = True# Flags indicating proximity to TC center

addvars.addaxesshape      = False# Array of variables based on 
                                 #  the major and minor axes from 
                                 #  eigenvalue/vectors
addvars.addaxesshape1mmhr = False# As above for >1mm/hr pixels
addvars.addaxesshape1mmhr = False# As above for >1mm/hr pixels
addvars.addaxesshape10mmhr= False# As above for >10mm/hr pixels
addvars.addasymmetry      = False# Asymmetry shape parameter 
                                 #  (Zick et al. 2016)
addvars.addasymmetryc     = False# As above for convective pixels
addvars.addfragmentation  = False# Fragmentation shape parameter 
                                 #  (Zick et al. 2016)
addvars.addfragmentationc = False# As above for convective pixels
addvars.adddispersion     = False# Dispersion shape parameter 
                                 #  (Zick et al. 2016)
addvars.adddispersionc    = False# As above for convective pixels

addvars.addperimeter      = False# Distance around perimeter of 
                                 #  shape (alpha-shape method)

# Inputs for specific variables
  
# Grid spacing in degrees lon, lat (only required for 
#  area/volrainrate/shape)
addvars.dx = 0.1
addvars.dy = 0.1

# Directory and filename of TC data (only required for TC 
#  information)
addvars.dataTCdir = "/uufs/chpc.utah.edu/common/home/varble-group2/IBTrACS/"
addvars.fileTCid  = "IBTrACS.ALL.v04r00.nc"

# Shape metric indexes
addvars.minshapesize       = 30  # Minimum number of pixels for axes
addvars.minshapegood       = 0.8 # Minimum goodness of fit for axes
addvars.minshapesize10mmhr = 15  # Minimum number of pixels for >10mm/hr axes
addvars.minshapegood10mmhr = 0.55 # Minimum goodness of fit for >10mm/hr axes

#==================================================================
# Namelist for add_ERA5_vars
#==================================================================
# Do not change
def addERA5vars():
  pass
#==================================================================

# Subset (for certain range of dates) 
addERA5vars.ssdat = True
addERA5vars.date1 = "20150801"
addERA5vars.date2 = "20150810"
addERA5vars.ssobs = False
addERA5vars.obid1 = "103231"
addERA5vars.obid2 = "103233"
addERA5vars.sslst = False
addERA5vars.lstfn = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/james/data/UG_file_lists/fn_ug_all_80_2015.txt"

# Number of processes for parrallelization
addERA5vars.serialorparallel = 2
addERA5vars.njobs = 20

# Type of area or mean
addERA5vars.addctarea = False # Area centered on TIPS
addERA5vars.addctmean = True  # Mean of centered area
addERA5vars.addinarea = False # All points in calculated inflow
addERA5vars.addinmean = False # Add mean of inflow region

# Rain check or no rain check
addERA5vars.addrainchk   = True
addERA5vars.addnorainchk = False

# Full, anomaly, and/or normalized anomaly
addERA5vars.addfull  = True
addERA5vars.addanom  = True
addERA5vars.addanomn = True

# Moisture Variables desired
addERA5vars.addTCWV = False # ERA5 Total Column Water Vapor
addERA5vars.addSH18 = False # ERA5 1000-850 hPa Specific Humidity
addERA5vars.addSH84 = False # ERA5 800-400 hPa Specific Humidity
addERA5vars.addMF18 = False # ERA5 1000-850 hPa moisture flux convergence
addERA5vars.addMF84 = False # ERA5 800-400 hPa moisture flux convergence
addERA5vars.addMA18 = False # ERA5 1000-850 hPa moisture advection
addERA5vars.addMA84 = False # ERA5 800-400 hPa moisture advection
addERA5vars.addMC18 = False # ERA5 1000-850 hPa moisture convergence
addERA5vars.addMC84 = False # ERA5 800-400 hPa moisture convergence

# Kinematic variables desired
addERA5vars.addMS18 = False # ERA5 1000-850 hPa shear magnitude
addERA5vars.addMS17 = False # ERA5 1000-700 hPa shear magnitude
addERA5vars.addDS17 = False # ERA5 1000-700 hPa shear direction
addERA5vars.addMS84 = False # ERA5 800-400 hPa shear magnitude
addERA5vars.addMS82 = False # ERA5 800-200 hPa shear magnitude
addERA5vars.addMS65 = False # ERA5 650-500 hPa shear magnitude
addERA5vars.addMS14 = False # ERA5 1000-400 hPa shear magnitude
addERA5vars.addCV18 = False # ERA5 1000-850 hPa convergence
addERA5vars.addCV31 = False # ERA5 300-100 hPa convergence
addERA5vars.addVE18 = False # ERA5 1000-850 hPa vertical motion
addERA5vars.addVE84 = False # ERA5 800-400 hPa vertical motion

# ERA5 domain variables
addERA5vars.hda           = 2.5 # Half data area in degrees
addERA5vars.hoursbefore   = np.arange(48,0,-3) # Hours before (descending)
addERA5vars.hoursafter    = addERA5vars.hoursbefore[::-1] # Hours after (ascending)
addERA5vars.avgmissfrac   = 0.8 # fraction of domain missing before average is not carried out

# Directory and filenames of ERA5 data
addERA5vars.fileTCRWid  = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/clouds/TCRW/ERA5.TCRW."
addERA5vars.fileTCWVid  = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/moisture/TCWV/ERA5.TCWV.anomaly_"
addERA5vars.fileSH18id  = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/moisture/SPHU/ERA5.SPHU_1000-850hPamean.anomaly_"
addERA5vars.fileSH84id  = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/moisture/SPHU/ERA5.SPHU_800-400hPamean.anomaly_"
addERA5vars.fileMF18id  = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/convergence/ERA5.MFC_1000-850hPamean.anomaly_"
addERA5vars.fileMF84id  = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/convergence/ERA5.MFC_800-400hPamean.anomaly_"
addERA5vars.fileMA18id  = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/convergence/ERA5.MADV_1000-850hPamean.anomaly_"
addERA5vars.fileMA84id  = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/convergence/ERA5.MADV_800-400hPamean.anomaly_"
addERA5vars.fileMC18id  = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/convergence/ERA5.MCONV_1000-850hPamean.anomaly_"
addERA5vars.fileMC84id  = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/convergence/ERA5.MCONV_800-400hPamean.anomaly_"
addERA5vars.fileMS18id = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/shear/ERA5.MSHR_1000-850hPa.anomaly_"
addERA5vars.fileMS17id = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/shear/ERA5.MSHR_1000-700hPa.anomaly_"
addERA5vars.fileDS17id = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/shear/ERA5.DSHR_1000-700hPa."
addERA5vars.fileMS84id = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/shear/ERA5.MSHR_800-400hPa.anomaly_"
addERA5vars.fileMS82id = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/shear/ERA5.MSHR_800-200hPa.anomaly_"
addERA5vars.fileMS65id = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/shear/ERA5.MSHR_650-500hPa.anomaly_"
addERA5vars.fileMS14id = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/shear/ERA5.MSHR_1000-400hPa.anomaly_"
addERA5vars.fileCV18id = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/convergence/ERA5.CONV_1000-850hPamean.anomaly_"
addERA5vars.fileCV31id = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/convergence/ERA5.CONV_300-100hPamean.anomaly_"
addERA5vars.fileVE18id = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/vertmotion/ERA5.VERT_1000-850hPamean.anomaly_"
addERA5vars.fileVE84id = "/uufs/chpc.utah.edu/common/home/u0816744/zips2/ERA5_derived/vertmotion/ERA5.VERT_800-400hPamean.anomaly_"

#==================================================================
# End namelist
#==================================================================

