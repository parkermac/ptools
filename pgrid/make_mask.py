# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 08:36:02 2016

@author: PM5

Code to create an initial mask for a grid.

"""

from importlib import reload
import gfun; reload(gfun)
gridname, dir0, pgdir, gdir = gfun.gstart()

import numpy as np
import shutil
import os

import zfun

#%% load grid file

if False:
    # interactive selection
    print('\n%s\n' % '** Choose file to edit **')
    fn_list_raw = os.listdir(gdir)
    fn_list = []
    for item in fn_list_raw:
        if item[-3:] == '.nc':
            fn_list.append(item)
    Nfn = len(fn_list)
    fn_dict = dict(zip(range(Nfn), fn_list))
    for nfn in range(Nfn):
        print(str(nfn) + ': ' + fn_list[nfn])
    my_nfn = int(input('-- Input number -- '))
    fn = fn_dict[my_nfn]
else:
    fn = 'grid_m00_r00_s00_x00.nc'

in_fn = gdir + fn

#%% create the new file name
fn_new = gfun.increment_filename(fn, tag='_m')
out_fn = gdir + fn_new

#%% get the grid from NetCDF
import netCDF4 as nc
ds = nc.Dataset(in_fn)
plon = ds.variables['lon_psi_ex'][:]
plat = ds.variables['lat_psi_ex'][:]
Z = ds.variables['z'][:]
mask_rho_orig = ds.variables['mask_rho'][:]
ds.close()
plon_vec = plon[0,:]
plat_vec = plat[:,0]

#%% coastline

cmat = gfun.get_coast()

#%% create a boolean mask array (True where masked = land)

# note that this is the opposite of the ROMS convention
# where mask_rho = 1. over water, and 0. over land

ZM = Z > 0

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
    ZM[jj0, ii0] = False

#%% Save the output file

# create the new mask
mask_rho = mask_rho_orig.copy()
mask_rho[ZM == True] = 0

if not np.all(mask_rho == mask_rho_orig):
    print('Creating ' + out_fn)
    try:
        os.remove(out_fn)
    except OSError:
        pass # assume error was because the file did not exist
    shutil.copyfile(in_fn, out_fn)
    ds = nc.Dataset(out_fn, 'a')
    ds['mask_rho'][:] = mask_rho
    ds.close()
else:
    print('No change to mask')
