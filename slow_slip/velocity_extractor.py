#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This creates and a single NetCDF file containing bottom velocity
for some time range.

"""

from datetime import datetime, timedelta
start_time = datetime.now()
import netCDF4 as nc

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import numpy as np
import zrfun
import zfun

pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

gridname = 'cascadia1'
tag = 'base'
ex_name = 'lobio5'
list_type = 'hourly'
#
Ldir = Lfun.Lstart(gridname, tag)
Ldir['ex_name'] = ex_name
Ldir['gtagex'] = Ldir['gtag'] + '_' + Ldir['ex_name']
#
if Ldir['lo_env'] == 'pm_mac':
    whichyear = 2017
    dt0 = datetime(whichyear,8,5)
    dt1 = datetime(whichyear,8,9)
elif Ldir['lo_env'] == 'pm_fjord':
    whichyear = 2017
    dt0 = datetime(whichyear,1,1)
    dt1 = datetime(whichyear,12,31)
#
if list_type == 'hourly':
    fn_list = []
    dt = dt0
    while dt <= dt1:
        date_string = dt.strftime(format='%Y.%m.%d')
        Ldir['date_string'] = date_string
        f_string = 'f' + Ldir['date_string']
        in_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/' + f_string + '/'
        if (list_type == 'low_passed') and ('low_passed.nc' in os.listdir(in_dir)):
            fn_list.append(in_dir + 'low_passed.nc') 
        elif list_type == 'hourly':
            for hh in range(2,26):
                hhhh = ('0000' + str(hh))[-4:]
                fn_list.append(in_dir + 'ocean_his_' + hhhh + '.nc')
        dt = dt + timedelta(days=1)
#
# make some things
fn = fn_list[0]
ds = nc.Dataset(fn)
G = zrfun.get_basic_info(fn, only_G=True)
h = ds['h'][:]
ds.close()
    
NT = len(fn_list)

itag = Ldir['gtagex'] + '_' + str(whichyear)

if Ldir['lo_env'] == 'pm_mac':
    etag = '_mac'
elif Ldir['lo_env'] == 'pm_fjord':
    etag = ''

# prepare a directory for results
outdir0 = Ldir['parent'] + 'ptools_output/slow_slip/'
Lfun.make_dir(outdir0, clean=False)
outdir = outdir0 + 'bottom_velocity_extractions/'
Lfun.make_dir(outdir, clean=False)
# output file
out_name = 'velocity_' + itag + etag + '.nc'
out_fn = outdir + out_name
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist
ds2 = nc.Dataset(out_fn, 'w')

# lists of variables to process
dlist = ['xi_rho', 'eta_rho', 'xi_psi', 'eta_psi', 'ocean_time']
vn_list2 = [ 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi', 'mask_rho', 'h']

ds1 = nc.Dataset(fn_list[0])

# make a blank for writing interpolated u and v onto
h = ds1['h'][:]
uvblank = 0. * h
mask = ds1['mask_rho'][:]
uvblank = np.ma.masked_where(mask==0, uvblank)

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
for vn in ['ocean_time']:
    varin = ds1[vn]
    vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
    vv.long_name = varin.long_name
    vv.units = varin.units
    try:
        vv.time = varin.time
    except AttributeError:
        # ocean_time has no time
        pass

ds1.close()

# add derived variables
vv = ds2.createVariable('u', float, ('ocean_time', 'eta_rho', 'xi_rho'))
vv.long_name = 'bottom eastward velocity'
vv.units = 'm s-1'
vv.time = 'ocean_time'
#
vv = ds2.createVariable('v', float, ('ocean_time', 'eta_rho', 'xi_rho'))
vv.long_name = 'bottom northward velocity'
vv.units = 'm s-1'
vv.time = 'ocean_time'


# copy data
NT = len(fn_list)
tt = 0
for fn in fn_list:
    ds = nc.Dataset(fn)
    ds2['ocean_time'][tt] = ds['ocean_time'][0].squeeze()
    # make u and v interpolated onto the rho grid
    u0 = ds['u'][0, 0, :, :].squeeze()
    v0 = ds['v'][0, 0, :, :].squeeze()
    u = uvblank.copy()
    v = uvblank.copy()
    u[:, 1:-1] = (u0[:, 1:] + u0[:, :-1])/2
    v[1:-1, :] = (v0[1:, :] + v0[:-1, :])/2
    ds2['u'][tt,:,:] = u
    ds2['v'][tt,:,:] = v

    tt += 1
    ds.close()

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

if Ldir['lo_env'] == 'pm_mac':
    # testing
    zfun.ncd(out_fn)
    #
    import matplotlib.pyplot as plt
    plt.close('all')
    fig = plt.figure()
    #
    ax = fig.add_subplot(121)
    ds = nc.Dataset(out_fn)
    vn = 'u'
    pch = ax.pcolormesh(ds['lon_psi'][:], ds['lat_psi'][:], ds[vn][-1,1:-1,1:-1].squeeze())
    fig.colorbar(pch)
    ax.set_title(ds[vn].long_name + ' (' + ds[vn].units + ')')
    pfun.dar(ax)
    #
    ax = fig.add_subplot(122)
    ds = nc.Dataset(out_fn)
    vn = 'v'
    pch = ax.pcolormesh(ds['lon_psi'][:], ds['lat_psi'][:], ds[vn][-1,1:-1,1:-1].squeeze())
    fig.colorbar(pch)
    ax.set_title(ds[vn].long_name + ' (' + ds[vn].units + ')')
    pfun.dar(ax)
    plt.show()


