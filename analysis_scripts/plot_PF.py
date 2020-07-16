#==========================================================
# Script fits an ellipse to IMERG precip rates for a given 
#  IPF and plots those precip rates along with the fitted 
#  ellipse.
# 
# Uses the least squares ellipse fitting method from:
#  https://github.com/bdhammel/least-squares-ellipse-fitting
# 
# James Russell 2019
#==========================================================

# Import libraries
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import nclcmaps
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Read data from current IPF netcdf file
f = Dataset("/uufs/chpc.utah.edu/common/home/varble-group2/james/FiT_CPEX-AW/TIPS_test/TIPS_100703_201808031200_10_-17.nc") 
datalat = f.groups["lats"].groups["data"]
datalon = f.groups["lons"].groups["data"]
datapcp = f.groups["instrain"].groups["data"]
datakeys = datalat.__dict__.keys()

# Loop over all times in IPF
for k in datakeys:
  print("Working on "+str(k))

  # Get data for current time
  if hasattr(datalon.getncattr(k), "__len__"):
    print("Plotting data")
    data = [datalon.getncattr(k),datalat.getncattr(k)]
    pcp  = datapcp.getncattr(k)

    # Make figure
    plt.close('all')
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.axis('equal')

    # Plot precipitation data
    sc = ax.scatter(data[0], data[1], c=pcp, marker=',',
         cmap=nclcmaps.cmap('WhiteBlueGreenYellowRed'),
         vmin=0.5,vmax=16.5)

    cb = plt.colorbar(sc, ax=ax)

    # Finalize and save figure
    plt.gca().set_aspect('equal', adjustable='box')
    cb.set_label('IMERG rain rate (mm/hr)')
    plt.xlim(np.mean(data[0])-3., np.mean(data[0])+3.)
    plt.ylim(np.mean(data[1])-3., np.mean(data[1])+3.)
    plt.annotate("IPF @ time: "+str(k), xy=(0.01, 0.97), xycoords='axes fraction')
    plt.savefig('figs/IMERGfit_'+str(k)+'.png')
    plt.close()


