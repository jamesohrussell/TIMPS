#==========================================================
# Script reads IPF data in a folder and plots all tracks on
#   a map of the tracked domain.
#
# James Russell 2019
#==========================================================

# Import python libraries (do not change)
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt; plt.rcdefaults()
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

# For plotting
pltxmn = 0   # Min lat for plotting
pltxmx = 24   # Max lat for plotting
pltymn = 0  # Min lon for plotting
pltymx = 2500  # Max lon for plotting

#==========================================================
# Read data
#==========================================================

# Open file
data = Dataset("ugdata_dts_2014-2018_ITCZ.nc")

# Read variables
bins = data.variables["bins"][:]
tims = data.variables["tims"][:]

#==========================================================
# Make plot
#==========================================================

print("Plotting data")

fig, ax = plt.subplots()
ax.bar(tims, bins, alpha=0.5, width=0.5,edgecolor='black')
ax.xaxis.set_major_locator(MultipleLocator(4))
ax.xaxis.set_minor_locator(MultipleLocator(1))
plt.xlabel('Local time')
plt.ylabel('# of occurences upscale growth')
plt.xlim(pltxmn, pltxmx)
plt.ylim(pltymn, pltymx)
plt.savefig('upgrdts.png')
