# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 16:28:41 2016

@author: PM5
"""

from importlib import reload
import gfun; reload(gfun)
G = gfun.gstart()

import netCDF4 as nc
import numpy as np
import shutil

import os

#%% select grid file
fn = gfun.select_file(G)
in_fn = G['gdir'] + fn
# create new file name
fn_new = gfun.increment_filename(fn, tag='_x')
out_fn = G['gdir'] + fn_new


#%% load the data

ds = nc.Dataset(in_fn)

mask_rho = ds.variables['mask_rho'][:]

ds.close()

#%% make the masks

mask_u_bool = (mask_rho[:, 1:] == 0) | (mask_rho[:, :-1] == 0)
mask_u = np.ones_like(mask_u_bool, dtype=int)
mask_u[mask_u_bool] = 0

mask_v_bool = (mask_rho[1:, :] == 0) | (mask_rho[:-1, :] == 0)
mask_v = np.ones_like(mask_v_bool, dtype=int)
mask_v[mask_v_bool] = 0

mask_psi_bool = ( (mask_rho[1:, 1:] == 0) | (mask_rho[:-1, :-1] == 0) |
                (mask_rho[1:, :-1] == 0) | (mask_rho[:-1, 1:] == 0) )
mask_psi = np.ones_like(mask_psi_bool, dtype=int)
mask_psi[mask_psi_bool] = 0

#%% save the output to NetCDF

# get rid of old version
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist
print('Creating ' + out_fn)
shutil.copyfile(in_fn, out_fn)

# open NetCDF file
ds = nc.Dataset(out_fn, 'a')

# add data to fields
tag_list = ['u', 'v', 'psi']
mask_dict = {'u': mask_u, 'v': mask_v, 'psi': mask_psi}
for tag in tag_list:
    ds['mask_'+tag][:] = mask_dict[tag]

ds.close()
