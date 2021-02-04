import xarray as xr
import numpy as np

ds = xr.open_dataset("IMERG_tracked_17614_4Dobjects.nc")

value =  np.squeeze(ds.value.values)

inds = np.nonzero(value)

print(set([int(v) for v in value[inds[0][:],inds[1][:]]]))
