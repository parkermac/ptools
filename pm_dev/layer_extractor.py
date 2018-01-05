#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 14:42:39 2017

@author: pm7

This creates and a single NetCDF file containing only selected fields
on a chosen sigma-surface, for some time range.

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
from datetime import datetime, timedelta
start_time = datetime.now()
import netCDF4 as nc
import numpy as np
import zrfun
import zfun

gridname = 'cascadia1'
tag = 'base'
ex_name = 'lobio1'
n_layer = 0
list_type = 'low_passed'
do_vel = False
out_name = 'ocean_layer.nc'

Ldir = Lfun.Lstart(gridname, tag)
Ldir['ex_name'] = ex_name
Ldir['gtagex'] = Ldir['gtag'] + '_' + Ldir['ex_name']

if Ldir['env'] == 'pm_mac':
    dt0 = datetime(2017,8,5)
    dt1 = datetime(2017,8,9)    
elif Ldir['env'] == 'pm_fjord':
    dt0 = datetime(2013,1,2)
    dt1 = datetime(2015,12,31)    

# prepare a directory for results
outdir0 = Ldir['parent'] + 'ptools_output/layer/'
Lfun.make_dir(outdir0, clean=False)
outdir = outdir0 + Ldir['gtagex'] + '/'
Lfun.make_dir(outdir, clean=True)

# output files
out_fn = outdir + out_name
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist

#%% make a list of files to extract from

if list_type == 'low_passed':    
    fn_list = []
    dt = dt0    
    while dt <= dt1:        
        date_string = dt.strftime(format='%Y.%m.%d')
        Ldir['date_string'] = date_string
        f_string = 'f' + Ldir['date_string']
        in_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string + '/'
        if 'low_passed.nc' in os.listdir(in_dir):
            fn_list.append(in_dir + 'low_passed.nc') 
        else:
            print('Missing file for ' + date_string)
        dt = dt + timedelta(days=1)

#%% create the layer NetCDF file

# make z
fn = fn_list[0]
ds = nc.Dataset(fn)
S = zrfun.get_basic_info(fn, only_S=True)
h = ds['h'][:]
z = zrfun.get_z(h, 0*h, S, only_rho=True)
z0 = z[n_layer,:,:].squeeze()
ds.close()

#%% Initialize the multi-file input dataset
# would this work for a whole year?
ds1 = nc.MFDataset(fn_list)

#%% make layer velocity
if do_vel:
    u0 = ds1['u'][:, n_layer, :, :].squeeze()
    v0 = ds1['v'][:, n_layer, :, :].squeeze()
    u = np.nan * ds1['salt'][:, n_layer, :, :].squeeze()
    v = u.copy()
    u[:, :, 1:-1] = (u0[:, :, 1:] + u0[:, :, :-1])/2
    v[:, 1:-1, :] = (v0[:, 1:, :] + v0[:, :-1, :])/2
    
ds2 = nc.Dataset(out_fn, 'w')

# lists of variables to process
dlist = ['xi_rho', 'eta_rho', 'xi_psi', 'eta_psi', 'ocean_time']
vn_list2 = [ 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi', 'mask_rho', 'h']
vn_list2t = ['zeta', 'ocean_time']
vn_list3t = ['salt', 'temp']
#vn_list2t = ['Uwind', 'Vwind', 'zeta', 'ocean_time']
#vn_list3t = ['salt', 'temp', 'NO3', 'phytoplankton',
#           'zooplankton', 'oxygen', 'TIC', 'alkalinity', 'PH', 'ARAG']

# Copy dimensions
for dname, the_dim in ds1.dimensions.items():
    if dname in dlist:
        ds2.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
    
# Copy variables
for vn in vn_list2:
    varin = ds1[vn]
    vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
    vv.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    vv[:] = ds1[vn][:]
#
for vn in vn_list2t:
    varin = ds1[vn]
    vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
    vv.long_name = varin.long_name
    vv.units = varin.units
    try:
        vv.time = varin.time
    except AttributeError:
        # ocean_time has no time
        pass
    vv[:] = ds1[vn][:]
#
for vn in vn_list3t:
    do_var = True
    try:
        varin = ds1[vn]
    except IndexError:
        # designed so that it still works when we don't have a variable from this list
        # e.g. when there is no bio or carbon
        do_var = False
        print(' - Variable not found: ' + vn)
    if do_var==True:
        dd = tuple([d for d in varin.dimensions if d != 's_rho'])
        vv = ds2.createVariable(vn, varin.dtype, dd)
        if vn == 'PH':
            vv.long_name = 'pH'
        elif vn == 'ARAG':
            vv.long_name = 'Aragonite Saturation State'
        else:
            vv.long_name = varin.long_name
        try:
            vv.units = varin.units
        except AttributeError:
            # salt has no units
            pass
        vv.time = varin.time
        vv[:] = ds1[vn][:, n_layer, :, :].squeeze()

# Add derived variables
vv = ds2.createVariable('z', float, ('eta_rho', 'xi_rho'))
vv.long_name = 'z position closest to free surface for 3D variables '
vv.units = 'meter'
vv[:] = z0
#
if do_vel:
    vv = ds2.createVariable('u', float, ('ocean_time', 'eta_rho', 'xi_rho'))
    vv.long_name = 'eastward near-surface velocity'
    vv.units = 'meter second-1'
    vv.time = 'ocean_time'
    vv[:] = u
    #
    vv = ds2.createVariable('v', float, ('ocean_time', 'eta_rho', 'xi_rho'))
    vv.long_name = 'northward near-surface velocity'
    vv.units = 'meter second-1'
    vv.time = 'ocean_time'
    vv[:] = v

ds1.close()
ds2.close()

#%% prepare for finale
import collections
result_dict = collections.OrderedDict()
time_format = '%Y.%m.%d %H:%M:%S'
result_dict['start_time'] = start_time.strftime(time_format)
end_time = datetime.now()
result_dict['end_time'] = end_time.strftime(time_format)
dt_sec = (end_time - start_time).seconds
result_dict['total_seconds'] = str(dt_sec)
if os.path.isfile(out_fn):
    result_dict['result'] = 'success'
else:
    result_dict['result'] = 'fail'
    
zfun.ncd(out_fn)


