# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 08:36:02 2016

@author: PM5

Code to create an initial mask for a grid.

"""

from importlib import reload
import gfun; reload(gfun)
Gr =gfun.gstart()
import pfun

import numpy as np
import shutil
import pickle
import os

import zfun

#% select grid file
fn = gfun.select_file(Gr)
in_fn = Gr['gdir'] + fn
# create new file name
fn_new = gfun.increment_filename(fn, tag='_m')
out_fn = Gr['gdir'] + fn_new

#% get the grid from NetCDF
import netCDF4 as nc
ds = nc.Dataset(in_fn)
plon = ds.variables['lon_psi_ex'][:]
plat = ds.variables['lat_psi_ex'][:]
z = -ds.variables['h'][:]
mask_rho_orig = ds.variables['mask_rho'][:]
ds.close()
plon_vec = plon[0,:]
plat_vec = plat[:,0]

# load the default choices
dch = pickle.load(open(Gr['gdir'] + 'choices.p', 'rb'))

# PROCESSING

# Create a boolean mask array (True where masked = land)
# following the numpy masked array convention.
# Note that this is the opposite of the ROMS convention
# where mask_rho = 1. over water, and 0. over land.
if mask_rho_orig.all() == 1:    
    print('Original mask all ones')
    # set z position of initial dividing line (positive up)
    m = z >= dch['z_land']
elif dch['do_cell_average']:
    # This branch applies when we created the grid using
    # do_cell_ave = True
    print('Original mask not all ones')
    #m = mask_rho_orig == 0 # seems OK (in low res)
    #m = mask_rho_orig < 1 # gets rid of a lot of water (in low res)
    #m = mask_rho_orig < .5 # OK but still removes a lot of PS
    m = mask_rho_orig < dch['z_land_alt'] # 0.1 maybe best?

# unmask the coast
if dch['unmask_coast']:
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
      
# remove islands and lakes
if dch['remove_islands']:
    # What this does is mask any water point that has land on 3 sides
    # or any land point that has water on three sides. By doing this repeatedly
    # you get rid of stray channels or peninsulas.
    # The number in range() determines how long of a feature is removed.
    # What the algorithm will not do, for example, is get rid of
    # a square lake of 4 cells.    
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
        
# make sure that any remaining high spots are masked
m[(m==False) & (z >= dch['z_land'])] = True

#% Save the output file

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
