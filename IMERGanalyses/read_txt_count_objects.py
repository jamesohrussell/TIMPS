#==========================================================
# Script to read FiT text output file and count number of 
# objects tracked
#
# James Russell 2019
#==========================================================

# Import python libraries
import numpy as np
import csv
from collections import Counter

#==========================================================
# Namelist
#==========================================================

# Directory and filename for FiT .txt ouput
datadir = "/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_100719/FiT_output_files_2017/"
fileid  = "IMERG_tracked_4Dobject_tree.txt"

# Min. hours object can last for
minhrs   = 6

# Observations per hour
obsperhr = 2

#==========================================================
# Function to read data
#==========================================================

def read_FiT_txt(dirandfile):
  "This function opens the FiT .txt output file and reads in the data"

  # Open file and give it the handle f
  f = open(dirandfile,'r') 

  # Use numpy function to read data from .eol file
  data = np.genfromtxt(f, dtype="str",delimiter="\t",skip_header=1)

  return(data)

  # Close file
  f.close()

#==========================================================
# Read data and count objects
#==========================================================

# Read data
data = read_FiT_txt(datadir+fileid)

# Create dictionary with number of occurrences of each IPF
obids = Counter([int(i) for i in data[:,0]])

# Count number of observations that last for 6 hours or more
obsgeminhrs = dict((k, v) for k, v in obids.items() if v >= minhrs*obsperhr)

# Print number of files 
print(str(len(obsgeminhrs))+" of "+str(len(obids))+" objects have duration >= "+str(minhrs)+" hours")

#==========================================================
# End Program
#==========================================================
