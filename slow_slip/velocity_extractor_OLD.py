#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This creates and a single NetCDF file containing bottom velocity,
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

n_layer = 0 # 0 = deepest
do_vel = False

gridname = 'cascadia1'
tag = 'base'
ex_name = 'lobio1'
list_type = 'low_passed'
#
Ldir = Lfun.Lstart(gridname, tag)
Ldir['ex_name'] = ex_name
Ldir['gtagex'] = Ldir['gtag'] + '_' + Ldir['ex_name']
tag = Ldir['gtagex']
#
if Ldir['env'] == 'pm_mac':
    dt0 = datetime(2017,8,5)
    dt1 = datetime(2017,8,9)
elif Ldir['env'] == 'fjord':
    dt0 = datetime(2013,1,2)
    dt1 = datetime(2015,12,31)
#
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
#
# make some things
fn = fn_list[0]
ds = nc.Dataset(fn)
G, S, T = zrfun.get_basic_info(fn)
h = ds['h'][:]
z = zrfun.get_z(h, 0*h, S, only_rho=True)
z0 = z[n_layer,:,:].squeeze()
ds.close()
    
if Ldir['env'] == 'pm_mac':
    etag = '_mac'
elif Ldir['env'] == 'fjord':
    etag = ''
    
NT = len(fn_list)

# prepare a directory for results
outdir0 = Ldir['parent'] + 'ptools_output/slow_slip/'
Lfun.make_dir(outdir0, clean=False)
outdir = outdir0 + 'bottom_velocity_extractions/'
Lfun.make_dir(outdir, clean=False)
# output file
out_name = 'velocity_' + tag + etag + '.nc'
out_fn = outdir + out_name
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist

#%% make bottom velocity

if do_vel:
    u = (np.nan*h) * np.ones((NT,1,1))
    v = u.copy()
    for tt in range(NT):
        fn = fn_list[tt]
        ds1 = nc.Dataset(fn)
        u0 = ds1['u'][0, 0, :, :].squeeze()
        v0 = ds1['v'][0, 0, :, :].squeeze()
        u = np.nan * ds1['salt'][0, 0, :, :].squeeze()
        v = u.copy()
        u[tt, :, 1:-1] = (u0[:, 1:] + u0[:, :-1])/2
        v[tt, 1:-1, :] = (v0[1:, :] + v0[:-1, :])/2
        ds1.close()

# initialize output Dataset
ds2 = nc.Dataset(out_fn, 'w')

# lists of variables to process
if model_type == 'LiveOcean':
    dlist = ['xi_rho', 'eta_rho', 'xi_psi', 'eta_psi', 'ocean_time']
    vn_list2 = [ 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi', 'mask_rho', 'h']
vn_list2t = ['zeta', 'ocean_time']
vn_list3t = ['salt', 'temp']

ds1 = nc.Dataset(fn_list[0])

# Copy dimensions
for dname, the_dim in ds1.dimensions.items():
    if dname in dlist:
        ds2.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

# Copy variables
for vn in vn_list2:
    if model_type == 'LiveOcean':
        varin = ds1[vn]
    vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
    vv.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    if model_type == 'LiveOcean':
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

# copy data
tt = 0
for fn in fn_list:
    ds = nc.Dataset(fn)
    ds2['ocean_time'][tt] = ds['ocean_time'][0].squeeze()
    ds2['zeta'][tt,:,:] = ds['zeta'][0, :, :].squeeze()
    tt += 1
    ds.close()
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
        vv.long_name = varin.long_name
        try:
            vv.units = varin.units
        except AttributeError:
            # salt has no units
            pass
        vv.time = varin.time
        for tt in range(NT):
            fn = fn_list[tt]
            ds = nc.Dataset(fn)
            vv[tt,:,:] = ds[vn][0, n_layer, :, :].squeeze()
            ds.close()

# Add derived variables
vv = ds2.createVariable('z', float, ('eta_rho', 'xi_rho'))
vv.long_name = 'z position closest to free surface for 3D variables '
vv.units = 'meter'
vv[:] = z0
#
if do_press:
    vv = ds2.createVariable('bp', float, ('ocean_time', 'eta_rho', 'xi_rho'))
    vv.long_name = 'bottom pressure'
    vv.units = 'Pa'
    vv.time = 'ocean_time'
    vv[:] = bp_arr
    #
    vv = ds2.createVariable('bpa', float, ('ocean_time', 'eta_rho', 'xi_rho'))
    vv.long_name = 'bottom pressure anomaly'
    vv.units = 'Pa'
    vv.time = 'ocean_time'
    vv[:] = bp_anom
    

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

if True:
    zfun.ncd(out_fn)
    
if Ldir['env'] == 'pm_mac':
    import matplotlib.pyplot as plt
    plt.close('all')
    
    fig = plt.figure()
    ax = fig.add_subplot(121)
    ds = nc.Dataset(out_fn)
    vn = 'bpa'
    pch = ax.pcolormesh(ds['lon_rho'][:], ds['lat_rho'][:], ds[vn][-1,:,:].squeeze())
    fig.colorbar(pch)
    ax.set_title(ds[vn].long_name + ' (' + ds[vn].units + ')')
    pfun.dar(ax)
    
    # put data into a DataFrame
    import pandas as pd
    zm = ds['zeta'][:,30,30].squeeze()
    tm = ds['ocean_time'][:]
    dtm_list = []
    for t in tm:
        dtm_list.append(Lfun.modtime_to_datetime(t))
    dti = pd.to_datetime(dtm_list)
    dti = dti.tz_localize('UTC')
    df = pd.DataFrame(data={'eta':zm}, index = dti)
    df.index.name = 'Date'
    ax = fig.add_subplot(122)
    df.plot(ax=ax)
    
    
    plt.show()
    ds.close()


