"""
Code to test isolayer plotting.
"""

# setup
import os
import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()
import zrfun
import zfun

plp = os.path.abspath('../../LiveOcean/plotting')
if plp not in sys.path:
    sys.path.append(plp)
import pfun

import numpy as np
import netCDF4 as nc

fn = '/Users/pm7/Documents/LiveOcean_roms/output/cas3_v0_lo6m/f2017.03.29/ocean_his_0002.nc'
ds = nc.Dataset(fn)

G, S, T = zrfun.get_basic_info(fn)

h = ds['h'][:]
zeta = 0*h
zw = zrfun.get_z(h, zeta, S, only_w=True)
DZ = np.diff(zw, axis=0)

mask = ds['mask_rho'][:]

s = ds['salt'][:].squeeze()
ss = 30

dz = DZ.copy()
dz[s > ss] = 0

zz = -dz.sum(axis=0)

zz[mask==0] = np.nan

zz[zz==0] = np.nan

zzm = np.ma.masked_where(np.isnan(zz), zz)

import matplotlib.pyplot as plt

#aa = [-123, -122.55, 47.5, 47.9]
aa = [-124, -122, 47, 49]

plt.close('all')
fig = plt.figure()
ax = fig.add_subplot(111)

cs = ax.pcolormesh(ds['lon_psi'][:], ds['lat_psi'][:], zzm[1:-1,1:-1],
                   vmin=-40, vmax=0, cmap='rainbow')
fig.colorbar(cs)
pfun.add_coast(ax)
ax.axis(aa)
pfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

plt.show()