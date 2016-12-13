# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 08:36:02 2016

@author: PM5

Code to create an initial mask for a grid.

"""

from importlib import reload
import gfun; reload(gfun)
G = gfun.gstart()
import pfun

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

#%% USER CHOICES

z_land = 0 # z position of initial dividing line (positive up)

unmask_coast = False

z_shallowest = -5
enforce_z_shallowest = False

remove_islands = True

#%% processing

# create a boolean mask array (True where masked = land)
m = z >= z_land
# note that this is the opposite of the ROMS convention
# where mask_rho = 1. over water, and 0. over land

# unmask the coast
if unmask_coast:
    # This unmasks it in the places where the
    # coastline crosses a tile, to facilitate wetting-drying
    cx, cy = pfun.get_coast()
    cmask = np.isnan(cx)
    cx = cx[~cmask]
    cy = cy[~cmask]
    ii0, ii1, ifr = zfun.get_interpolant(cx, plon_vec, extrap_nan=True)
    jj0, jj1, jfr = zfun.get_interpolant(cy, plat_vec, extrap_nan=True)
    # Don't unmask extrapolated points.
    ii0 = ii0[~np.isnan(ifr) & ~np.isnan(jfr)]
    jj0 = jj0[~np.isnan(ifr) & ~np.isnan(jfr)]
    m[jj0, ii0] = False

# enforce a minimum depth in unmasked cells
# (don't use with wet-dry applications)
if enforce_z_shallowest:
    z[(~m) & (z>z_shallowest)] = z_shallowest
      
# remove islands
if remove_islands:
    for ii in range(5):
        NR, NC = m.shape
        mm = m[1:-1, 1:-1]
        mn = m[2:, 1:-1]
        ms = m[:-2, 1:-1]
        me = m[1:-1, 2:]
        mw = m[1:-1, :-2]
        # remove islands of ocean
        MMo = ~mm & mn & ms & me
        mm[MMo] = True
        MMo = ~mm & mn & ms & mw
        mm[MMo] = True
        MMo = ~mm & mn & me & mw
        mm[MMo] = True
        MMo = ~mm & ms & me & mw
        mm[MMo] = True
        # remove islands of land
        MMl = mm & ~mn & ~ms & ~me
        mm[MMl] = False
        MMl = mm & ~mn & ~ms & ~mw
        mm[MMl] = False
        MMl = mm & ~mn & ~me & ~mw
        mm[MMl] = False
        MMl = mm & ~ms & ~me & ~mw
        mm[MMl] = False
        m[1:-1, 1:-1] = mm
       

#%% Save the output file

# create the new mask_rho
# 1 = water
# 0 = land
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
