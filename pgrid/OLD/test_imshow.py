"""
Test of imshow for edit_mask.py.

"""

import gfun
from importlib import reload
reload(gfun)
G = gfun.gstart()
# running gfun.gstart() sets the path to include pfun and zfun
import pfun
import zfun

import netCDF4 as nc

import numpy as np
import shutil
import os
import pandas as pd
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import matplotlib.colors as pltc
from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

gname = 'salish1/grid_m02_r01_s00_x00.nc'
in_fn = '/Users/PM5/Documents/ptools_output/pgrid/' + gname

# Test: retrieve the output

ds = nc.Dataset(in_fn)

H = ds.variables['h'][:]
mask_rho = ds.variables['mask_rho'][:]
plon = ds.variables['lon_psi_ex'][:]
plat = ds.variables['lat_psi_ex'][:]
lon = ds.variables['lon_rho'][:]
lat = ds.variables['lat_rho'][:]

ds.close()

# flip to work with imshow
h = np.flipud(H)
m = np.flipud(mask_rho)
# mask_rho:
# 1 = water
# 0 = land
hh = h.copy()
hh[m==0] = np.nan

# plotting

plt.close()

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111)

a = ax.matshow(h, interpolation='none', vmin=-100, vmax = 100)

plt.show()

for ii in range(5):
    
    x = plt.ginput(n=1)
    xx = np.array(x, dtype=int)
    if xx.shape != (1,2):
        xx = np.array([[-1000, -1000]])
    ix = xx[0, 0]
    iy = xx[0, 1]
    
    if np.isnan(hh[iy, ix]):
        hh[iy, ix] = h[iy, ix]
    else:
        hh[iy, ix] = np.nan

    a.set_data(hh)
    
    plt.draw()




