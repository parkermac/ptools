"""
Testing the NCO NetCDF operators for processing of a surface extraction.

RESULT: ncwa is actually slower in this case than doing it by hand, but
even more troubling I cannot see any difference between the plain average
from ncra and the "weighted" average using ncwa, so I really don't
understand what is happening.

Note the non-pythonic numbering used in the -d hyperslab notation
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from time import time
import subprocess

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
import zfun

sys.path.append(os.path.abspath('../../LiveOcean/plotting'))
import pfun

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

Ldir = Lfun.Lstart()
in_dir = Ldir['LOo'] + 'layer/cas6_v3_lo8b_2019.06.01_2019.08.31/'
in_name = 'surface_hourly_weight.nc'
in_fn = in_dir + in_name
# a total of 92 days, so we can form 90 tidal averages I believe

gw = zfun.godin_shape()
LW = len(gw)
ds = nc.Dataset(in_fn,'a')
if 'godin_weights' not in ds.variables:
    ds.createDimension('gw_dim', LW)
    vv = ds.createVariable('godin_weights', float, ('gw_dim'))
    vv.long_name = 'Godin Filter Weights'
    vv[:] = gw
else:
    # assume godin_weights already exists
    pass
ds.close()

# make a tidal average using ncwa
out_fn = in_dir + 'nco_test.nc'
vstr = 'zeta,lon_psi,lat_psi'
if False:
    cmd_list = ['ncwa',
        '-a', 'ocean_time',
        '-d', 'ocean_time,0,70',
        '-v',vstr,
        '-w','godin_weights',
        '-O', in_fn, out_fn]
else:
    cmd_list = ['ncra',
        '-d', 'ocean_time,0,70',
        '-v',vstr,
        '-O', in_fn, out_fn]
tt0 = time()
ret1 = subprocess.call(cmd_list)
print('Time to make average using NCO = %0.1f sec' % (time() - tt0))

ds = nc.Dataset(out_fn)
x = ds['lon_psi'][:]
y = ds['lat_psi'][:]
z = ds['zeta'][:].squeeze()
ds.close()

# make a tidal average by hand
tt0 = time()
zz = np.zeros_like(z)
ds = nc.Dataset(in_fn)
for ii in range(LW):
    #zz += gw[ii]*ds['zeta'][ii,:,:]
    zz += ds['zeta'][ii,:,:]/LW
ds.close()
print('Time to make average by hand = %0.1f sec' % (time() - tt0))

# PLOTTING
plt.close('all')
fig = plt.figure(figsize=(20,8))

ax = fig.add_subplot(131)
cs = ax.pcolormesh(x, y, z[1:-1,1:-1], vmin=-.1, vmax = .1, cmap='jet')
fig.colorbar(cs)
pfun.dar(ax)

ax = fig.add_subplot(132)
cs = ax.pcolormesh(x, y, zz[1:-1,1:-1], vmin=-.1, vmax = .1, cmap='jet')
fig.colorbar(cs)
pfun.dar(ax)

ax = fig.add_subplot(133)
cs = ax.pcolormesh(x, y, zz[1:-1,1:-1] - z[1:-1,1:-1], vmin=-.01, vmax = .01, cmap='jet')
fig.colorbar(cs)
pfun.dar(ax)

plt.show()



