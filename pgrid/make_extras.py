# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 16:28:41 2016

@author: PM5
"""

from importlib import reload
import gfun; reload(gfun)
G = gfun.gstart()

import netCDF4 as nc

import os
import seawater as sw

#%% select grid file
fn = gfun.select_file(G)
in_fn = G['gdir'] + fn
# create new file name
fn_new = gfun.increment_filename(fn, tag='_x')
out_fn = G['gdir'] + fn_new


#%% load the data

ds = nc.Dataset(in_fn)

lat = ds.variables['lat_rho'][:]
mask_rho = ds.variables['mask_rho'][:]

ds.close()

#%% make the masks

mask_u = mask_rho[:, 1:]
mask_v = mask_rho[1:, :]
mask_psi = mask_rho[1:, 1:]

#%% save the output to NetCDF

# get rid of old version
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist

# create new NetCDF file
foo = nc.Dataset(out_fn, 'w')

# create dimensions
M, L = lat.shape # use ROMS teminology
size_dict = {'rho': (M, L),
             'u': (M, L-1),
             'v': (M-1, L),
             'psi': (M-1, L-1),
             'psi_ex': (M+1, L+1)}
tag_list = ['rho', 'u', 'v', 'psi', 'psi_ex']
for tag in tag_list:
    foo.createDimension('eta_'+tag, size_dict[tag][0])
    foo.createDimension('xi_'+tag, size_dict[tag][1])

# create variables
tag_list = ['u', 'v', 'psi']
mask_var = dict()
for tag in tag_list:
    mask_var[tag] = foo.createVariable('lat_'+tag, float, ('eta_'+tag, 'xi_'+tag))
f_var = foo.createVariable('f', float, ('eta_rho', 'xi_rho'))


# add data to fields
mask_dict = {'u': mask_u, 'v': mask_v, 'psi': mask_psi}
for tag in tag_list:
    mask_var[tag][:] = mask_dict[tag]

f_var[:] = sw.f(lat)

foo.close()
