# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 08:36:02 2016

@author: PM5

Code to create an initial mask for a grid.

"""

from importlib import reload
import gfun; reload(gfun)
G = gfun.gstart()

import numpy as np
import shutil
import os

import zfun

#%% select grid file
fn = gfun.select_file(G)
in_fn = G['gdir'] + fn
# create new file name
fn_new = gfun.increment_filename(fn, tag='_m')
out_fn = G['gdir'] + fn_new

#%% get the grid from NetCDF
import netCDF4 as nc
ds = nc.Dataset(in_fn)
plon = ds.variables['lon_psi_ex'][:]
plat = ds.variables['lat_psi_ex'][:]
z = -ds.variables['h'][:]
mask_rho_orig = ds.variables['mask_rho'][:]
ds.close()
plon_vec = plon[0,:]
plat_vec = plat[:,0]

#%% coastline

cmat = gfun.get_coast()

#%% create a boolean mask array (True where masked = land)

# note that this is the opposite of the ROMS convention
# where mask_rho = 1. over water, and 0. over land

m = z > 0

if True:
    # This unmasks it in the places where the
    # coastline crosses a tile, to facilitate wetting-drying
    cx = cmat['lon']
    cy = cmat['lat']
    cmask = np.isnan(cx)
    cx = cx[~cmask]
    cy = cy[~cmask]
    ii0, ii1, ifr = zfun.get_interpolant(cx, plon_vec, extrap_nan=True)
    jj0, jj1, jfr = zfun.get_interpolant(cy, plat_vec, extrap_nan=True)
    # Don't unmask extrapolated points.
    ii0 = ii0[~np.isnan(ifr) & ~np.isnan(jfr)]
    jj0 = jj0[~np.isnan(ifr) & ~np.isnan(jfr)]
    m[jj0, ii0] = False

#%% enforce a minimum depth in unmasked cells
# we will need to change this for truly wet-dry applications
z[(~m) & (z>-5)] = -5

#%% Save the output file

# create the new mask_rho
mask_rho = np.ones_like(mask_rho_orig)
mask_rho[m == True] = 0

if not np.all(mask_rho == mask_rho_orig):
    print('Creating ' + out_fn)
    try:
        os.remove(out_fn)
    except OSError:
        pass # assume error was because the file did not exist
    shutil.copyfile(in_fn, out_fn)
    ds = nc.Dataset(out_fn, 'a')
    ds['mask_rho'][:] = mask_rho
    ds['h'][:] = -z
    ds.close()
else:
    print('No change to mask')
