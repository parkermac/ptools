# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 16:28:41 2016

@author: PM5
"""

from importlib import reload
import gfun; reload(gfun)
Gr =gfun.gstart()

import netCDF4 as nc
import numpy as np
import shutil

import os
import pickle

# load the default choices
dch = pickle.load(open(Gr['gdir'] + 'choices.p', 'rb'))

#%% select grid file
fn = gfun.select_file(Gr)
in_fn = Gr['gdir'] + fn
# create new file name
fn_new = gfun.increment_filename(fn, tag='_x')
out_fn = Gr['gdir'] + fn_new


#%% load the data

ds = nc.Dataset(in_fn, 'a')

mask_rho = ds['mask_rho'][:]

if dch['use_min_depth']:
    # Enforce a minimum depth (the default is that we ALWAYS do this).
    # And in this case it applies to the whole grid, not just the
    # masked part.  Probably not necessary, but it doesn't hurt, and makes
    # some plotting choices easier.  I wonder if it would make sense to
    # just do it at the very start...?
    h = ds['h'][:]
    hnew = h.copy()
    hnew[ h <= dch['min_depth'] ] = dch['min_depth']
    ds['h'][:] = hnew

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
    
# add global attributes (needed for matlab nesting code?)
ds.type = 'GRID file'
ds.history = 'whatever'
#
# add some grid info (needed for matlab nesting code?)
x_var = ds.createVariable('x_rho', float, ('eta_rho', 'xi_rho'))
y_var = ds.createVariable('y_rho', float, ('eta_rho', 'xi_rho'))
pm = ds['pm'][:]
pn = ds['pn'][:]
dx = 1/pm
dy = 1/pn
x_rho = np.cumsum(dx, axis=1)
y_rho = np.cumsum(dy, axis=0)
x_var[:] = x_rho
y_var[:] = y_rho
#
# and more
dndx_var = ds.createVariable('dndx', float, ('eta_rho', 'xi_rho'))
dmde_var = ds.createVariable('dmde', float, ('eta_rho', 'xi_rho'))
# Lp, Mp = x_rho.shape
# L = Lp-1
# Lm = L-1
# M = Mp-1
# Mm = M-1
dndx = np.zeros_like(x_rho)
dmde = np.zeros_like(x_rho)
# dndx(2:L,2:M) = 0.5.*(1.0./pn(3:Lp,2:M ) - 1.0./pn(1:Lm,2:M ));
# dmde(2:L,2:M) = 0.5.*(1.0./pm(2:L ,3:Mp) - 1.0./pm(2:L ,1:Mm));
dndx[1:-1,1:-1] = 0.5*(1.0/pn[2:,1:-1] - 1.0/pn[:-2,1:-1] )
dmde[1:-1,1:-1] = 0.5*(1.0/pm[1:-1,2:] - 1.0/pm[1:-1,:-2])
dndx_var[:] = dndx
dmde_var[:] = dmde

ds.close()
