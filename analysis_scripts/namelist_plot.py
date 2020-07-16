#==========================================================
# Namelist
#==========================================================

# Directory and filename for IMERG data (specify format --- only hdf5 or netcdf4 files)
# Should obtain from create_FiT_input_files.py
datadirIM = "/uufs/chpc.utah.edu/common/home/varble-group2/IMERG/"
fileidIM  = "3B-HHR.MS.MRG.3IMERG."
datahdf5  = True
datanc4   = False

# Date and time range
starttime = "20180101"
endtime   = "20180102" # Actual last day is day before

# Directory and filename for FiT output data files
datadirFi = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_CPEX-AW/FiT_output_test/"
fileidFi1 = "IMERG_tracked_"
fileidFi2 = "_4Dobjects.nc"
fileidFi3 = "IMERG_tracked_4Dobject_tree.txt"

# Directory and filename for IPFs
datadirPF = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_CPEX-AW/TIPS_test/"
fileidPF  = "TIPS_"

# Subset regions (if global set ssreg=False)
# Should obtain from create_FiT_input_files.py
ssreg  = True
latN   = 30
latS   = -10
lonW   = -70
lonE   = 0
ssname = "Atl"

# Gaussian smoothing (No smoothing if False)
# Should obtain from create_FiT_input_files.py
smooth = False
gaussian = False
stddev   = 1.5  # Standard deviation for gaussian smoothing
uniform  = False
width    = 5    # Number of points to spread running average

# if plotobj=True, plots object outlines 
# if idobj=False, outlines all objects after obid1
# if idobj=True, outlines only object with obid1
# if longobj=True, plots nlong longest lived objects
plotobj = False
idobj   = False
obid1   = 100039
longobj = False
nlong   = 10

# If plotthr=True plots thresholds used for FiT
# This option should be False if plotobj = True as this 
#   will just overwrite the contours there. 
# Should obtain thresholds from create_FiT_input_files.py
plotthr = False
tholds  = [0.8,3.0,9.0,27.0]

# For plotting
plot  = True # Turn on or off plot (mostly for debugging)
labnx = 4    # Number of x-labels between/incl. lonW,lonE
labny = 6    # Number of y-labels between/incl. latS,latN
dx    = 0.1  # Grid-spacing in x (0.1 for IMERG)
dy    = 0.1  # Grid spacing in y (0.1 for IMERG)
latmn = -20   # Min lat for plotting
latmx = 20   # Max lat for plotting
lonmn = -66 # Min lon for plotting
lonmx = 20  # Max lon for plotting
mincl = -0.5 # Min for colorbar
maxcl = 16.5 # Max for colorbar
#tracks= True # Plot tracks?

# For parallelization
njobs = 32
