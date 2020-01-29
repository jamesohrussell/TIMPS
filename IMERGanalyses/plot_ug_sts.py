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
import datetime as dt
from matplotlib import dates
from matplotlib.pyplot import figure

# For plotting
pltxmn = 0   # Min lat for plotting
pltxmx = 24   # Max lat for plotting
pltymn = 0.0  # Min lon for plotting
pltymx = 1.05  # Max lon for plotting

#==========================================================
# Read data
#==========================================================

# Open file
data = Dataset("ugdata_sts_2018_0-20N_30W-20W.nc")

# Read variables
bins = data.variables["bins"][:]
tims = data.variables["tims"][:]

dts = np.zeros((len(tims),6),dtype=int)
for i in range(len(tims)):
 dti = dt.datetime(2018, 1, 1) + dt.timedelta(tims[i] - 1)
 dts[i,0] = dti.year
 dts[i,1] = dti.month
 dts[i,2] = dti.day
 dts[i,3] = 0
 dts[i,4] = 0

dt1 = []
for i in range(len(dts)):
  dt1.append(str(dts[i,0])+"-"+str(dts[i,1])+"-"+str(dts[i,2])+" "+str(dts[i,3])+":"+str(dts[i,4]))
x = [int(i) for i in dates.datestr2num(dt1)]
formatter = dates.DateFormatter('%Y-%m-%d')

print(x)

#==========================================================
# Make plot
#==========================================================

print("Plotting data")

fig, ax = plt.subplots()
ax.bar(x, bins/max(bins), alpha=0.5, width=1.0,edgecolor='black')
plt.ylabel('Relative occurrence of upscale growth')
#plt.title('Upscale growth for the 2018 JJAS period')
plt.xticks(np.arange(min(x), max(x), 14))
ax = plt.gcf().axes[0]
ax.xaxis.set_major_formatter(formatter)
plt.gcf().autofmt_xdate(rotation=25)
plt.xlim(min(x), max(x))
ax.xaxis.set_major_locator(MultipleLocator(14))
ax.xaxis.set_minor_locator(MultipleLocator(7))
plt.ylim(pltymn, pltymx)
plt.savefig('upgrsts_2018.png')
