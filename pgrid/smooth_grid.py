# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 16:45:54 2016

@author: PM5

Smooth the grid.

"""

from importlib import reload
import gfun; reload(gfun)
G = gfun.gstart()
import pfun

import netCDF4 as nc
import os
import shutil
import numpy as np

#%% select grid file
fn = gfun.select_file(G)
in_fn = G['gdir'] + fn
# create new file name
fn_new = gfun.increment_filename(fn, tag='_s')
out_fn = G['gdir'] + fn_new

#%% Test: retrieve the output

ds = nc.Dataset(in_fn)
z = -ds.variables['h'][:]
mask_rho = ds.variables['mask_rho'][:]
plon = ds.variables['lon_psi_ex'][:]
plat = ds.variables['lat_psi_ex'][:]
lon = ds.variables['lon_rho'][:]
lat = ds.variables['lat_rho'][:]
dx = 1/ds.variables['pm'][:]
dy = 1/ds.variables['pn'][:]
ds.close()

#%% Smoothing set up

MSK = np.ones_like(lon).astype('float')
MSK[z>=0] = 0.
Hobs = -z.copy()
rx0max = 0.15

# make sure that anything not masked is
# not shallower than min_depth
min_depth = 0. # a negative number would be like on land
Hobs_1 = -z.copy()
Hobs = Hobs_1.copy()
Hobs_1[Hobs_1 < min_depth] = min_depth
Hobs[mask_rho==1] = Hobs_1[mask_rho==1]
# and make it positive definite under water
offset = 1
Hobs = Hobs + offset
# offset is removed later

#%% create the area matrix

AreaMatrix = dx * dy

#%% create smoothed bathymetry

import time
tt0 = time.time()
Hnew = gfun.GRID_PlusMinusScheme_rx0(MSK, Hobs, rx0max, AreaMatrix)
Hnew = Hnew - offset
print('Smoothing took %0.1f seconds' % (time.time() - tt0))
zn = -Hnew

# enforce a minimum depth
if False:
    zn[zn > -5.] = -5.

# Save the output file

print('Creating ' + out_fn)
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist
shutil.copyfile(in_fn, out_fn)
ds = nc.Dataset(out_fn, 'a')
ds['h'][:] = -zn

#%% plotting

if True:

    import matplotlib.pyplot as plt

    ax_lims = (plon[0,0], plon[0,-1], plat[0,0], plat[-1,0])

    z_m = np.ma.masked_where(mask_rho==0, z)
    zn_m = np.ma.masked_where(mask_rho==0, zn)
    dz = zn - z

    plt.close()

    fig = plt.figure(figsize=(14,10))

    ax = fig.add_subplot(121)
    cmap1 = plt.get_cmap(name='terrain')
    cs = ax.pcolormesh(plon, plat, zn_m,
                       vmin=-200, vmax=200, cmap = cmap1)
    fig.colorbar(cs, ax=ax, extend='both')
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(pfun.get_aa(ds))
    ax.set_title(G['gridname'] + '/' + fn_new)

    ax = fig.add_subplot(122)
    cmap1 = plt.get_cmap(name='bwr')
    cs = ax.pcolormesh(plon, plat, dz,
                       vmin=-100, vmax=100, cmap = cmap1)
    fig.colorbar(cs, ax=ax, extend='both')
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(pfun.get_aa(ds))
    ax.set_title('New Z - Old Z (m)')

    plt.show()

ds.close()
