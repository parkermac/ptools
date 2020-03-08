"""
Code to make a supercool plot of LiveOcean fields.
"""

import sys, os
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart()

sys.path.append(os.path.abspath('../../LiveOcean/plotting'))
import pfun

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

fn = Ldir['roms'] + 'output/cas6_v3_lo8b/f2019.07.04/ocean_his_0020.nc'
ds = nc.Dataset(fn)

x = ds['lon_rho'][:]
y = ds['lat_rho'][:]
h = ds['h'][:]
m = ds['mask_rho'][:]
z = -h
DO = ds['oxygen'][0,0,:,:]
DO = DO.data

X = x[::10,::10]
Y = y[::10,::10]
Z = z[::10,::10]

# PLOTTING
plt.close('all')
fig = plt.figure(figsize=(12,8))
ax = Axes3D(fig)

ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0, antialiased=False)#, facecolors=cm.jet(DO/10))
    
plt.show()


