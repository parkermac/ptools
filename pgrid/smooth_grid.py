# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 16:45:54 2016

@author: PM5

Smooth the grid.

"""

from importlib import reload
import gfun; reload(gfun)
G = gfun.gstart()

import netCDF4 as nc
import os
import shutil
import zfun
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
dx = ds.variables['dx'][:]
dy = ds.variables['dy'][:]
ds.close()

#%% Smoothing set up

MSK = np.ones_like(lon).astype('float')
MSK[z>=0] = 0.
Hobs = -z.copy()
rx0max = 0.15

#%% create the area matrix
#ulon = plon[:-1, :] + np.diff(plon, axis=0)/2
#vlat = plat[:, :-1]
#R = zfun.earth_rad(np.mean(plat[:,1]))
#dx = R * np.cos(np.pi*lat/180) * (np.pi * np.diff(ulon, axis = 1) / 180)
#dy = R * (np.pi * np.diff(vlat, axis = 0) / 180)
AreaMatrix = dx * dy

#%% create smoothed bathymetry
import time
tt0 = time.time()

if False:
    # slow version
    Hnew, info1, info2 = gfun.GRID_PlusMinusScheme_rx0(
                            MSK, Hobs, rx0max, AreaMatrix)
else:
    # fast version
    Hnew = gfun.GRID_PlusMinusScheme_rx0_v2(
                            MSK, Hobs, rx0max, AreaMatrix)

print('Smoothing took %0.1f seconds' % (time.time() - tt0))

zn = -Hnew

#%% Save the output file

print('Creating ' + out_fn)
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist
shutil.copyfile(in_fn, out_fn)
ds = nc.Dataset(out_fn, 'a')
ds['h'][:] = -zn
ds.close()

#%% plotting

if True:

    import matplotlib.pyplot as plt

    cmat = gfun.get_coast()

    ax_lims = (plon[0,0], plon[0,-1], plat[0,0], plat[-1,0])

    z_m = np.ma.masked_where(mask_rho==0, z)
    zn_m = np.ma.masked_where(mask_rho==0, zn)
    dz = zn - z

    plt.close()

    fig = plt.figure(figsize=(18,12))

    ax = fig.add_subplot(121)
    cmap1 = plt.get_cmap(name='terrain')
    cs = ax.pcolormesh(plon, plat, zn_m,
                       vmin=-200, vmax=200, cmap = cmap1)
    ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
    zfun.dar(ax)
    fig.colorbar(cs, ax=ax, extend='both')
    ax.set_xlim(ax_lims[:2])
    ax.set_ylim(ax_lims[-2:])
    ax.set_title(G['gridname'] + '/' + fn_new)

    ax = fig.add_subplot(122)
    cmap1 = plt.get_cmap(name='bwr')
    cs = ax.pcolormesh(plon, plat, dz,
                       vmin=-100, vmax=100, cmap = cmap1)
    ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
    zfun.dar(ax)
    fig.colorbar(cs, ax=ax, extend='both')
    ax.set_xlim(ax_lims[:2])
    ax.set_ylim(ax_lims[-2:])
    ax.set_title('New Z - Old Z (m)')

    plt.show()