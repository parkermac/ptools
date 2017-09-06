#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This creates and a single NetCDF file containing bottom pressure
and other fields, for some time range.

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
do_press = True

#model_type = 'LiveOcean'
model_type = 'Kurapov'

if model_type == 'LiveOcean':
    gridname = 'cascadia1'
    tag = 'base'
    ex_name = 'lobio1'
    list_type = 'low_passed'
    #
    Ldir = Lfun.Lstart(gridname, tag)
    Ldir['ex_name'] = ex_name
    Ldir['gtagex'] = Ldir['gtag'] + '_' + Ldir['ex_name']
    folder_tag = Ldir['gtagex']
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
    
elif model_type == 'Kurapov':
    Ldir = Lfun.Lstart()
    in_dir = Ldir['parent'] + 'ptools_data/Kurapov/'
    folder_tag = 'kurapov'
    
    # get grid info
    fng = in_dir + 'grd_wcofs_large_visc200.nc'
    dsg = nc.Dataset(fng)
    if False:
        print('\nGRID INFO')
        zfun.ncd(dsg)
    lon = dsg['lon_rho'][1030-1:1521, 375-1:615]
    lat = dsg['lat_rho'][1030-1:1521, 375-1:615]
    dsg.close()
        
    # make file list
    fn_list = []
    if Ldir['env'] == 'pm_mac':
        frange = range(1,5+1)
    elif Ldir['env'] == 'fjord':
        # we have 1-2276, so use range(1,2276+1)
        frange = range(1,2276+1)
    for ii in frange: 
        nn = ('0000' + str(ii))[-4:]
        fn_list.append(in_dir + 'Exp_29_Files/zts_ORWA_Parker_Exp29_' + nn + '.nc')
    
    # make some things
    S_info_dict = {'VTRANSFORM':2, 'VSTRETCHING':4, 'THETA_S':8,
        'THETA_B':3, 'TCLINE':50, 'N':40}
    S = zrfun.get_S(S_info_dict)
    S['s_rho'] = S['sc_r']
    S['s_w'] = S['sc_w']
    G = zrfun.get_basic_info(fng, only_G=True)
    fn = fn_list[0]
    ds = nc.Dataset(fn)
    h = ds['h'][:]
    z = zrfun.get_z(h, 0*h, S, only_rho=True)
    z0 = z[n_layer,:,:].squeeze()
    ds.close()
    
NT = len(fn_list)

# prepare a directory for results
outdir0 = Ldir['parent'] + 'ptools_output/slow_slip/'
Lfun.make_dir(outdir0, clean=False)
outdir = outdir0 + folder_tag + '/'
Lfun.make_dir(outdir, clean=True)
# output file
out_name = 'bottom_pressure.nc'
out_fn = outdir + out_name

# make bottom pressure
if do_press:
    g = 9.8
    for tt in range(NT):
        fn = fn_list[tt]
        ds1 = nc.Dataset(fn)
        if np.mod(tt,10)==0:
            print('tt = ' + str(tt) + '/' + str(NT) + ' ' + str(datetime.now()))
        if model_type == 'LiveOcean':
            zeta = ds1['zeta'][0,:,:].squeeze()
            if tt == 0:
                bp_arr = (0*zeta) * np.ones((NT,1,1))
                DA = G['DX'] * G['DY']
                DAm = np.ma.masked_where(zeta.mask, DA)
            rho = ds1['rho'][0,:,:,:].squeeze() + 1000.
            # note that rho appears to be in situ density, not potential density 
            z_w = zrfun.get_z(G['h'], zeta, S, only_w=True)
            DZ = np.diff(z_w, axis=0)
            bp_arr[tt,:,:] = (g * rho * DZ).sum(axis=0)
        elif model_type == 'Kurapov':
            import seawater
            zeta = ds1['zeta'][0,:,:].squeeze()
            if tt == 0:
                bp_arr = (0*zeta) * np.ones((NT,1,1))
                DA = G['DX'] * G['DY']
                DA = DA[1030-1:1521, 375-1:615]
                DAm = np.ma.masked_where(zeta.mask, DA)
            salt = ds1['salt'][0,:,:,:].squeeze()
            ptemp = ds1['temp'][0,:,:,:].squeeze() # potential temperature
            z_r = zrfun.get_z(G['h'][1030-1:1521, 375-1:615], 0* zeta, S, only_rho=True)
            z_w = zrfun.get_z(G['h'][1030-1:1521, 375-1:615], zeta, S, only_w=True)
            p = seawater.pres(-z_r, G['lat_rho'][0,0])
            temp = seawater.ptmp(salt, ptemp, 0*p, pr=p) # in situ temperature
            # for some reason seawater.dens throws errors if we don't do this
            sd = salt.data
            td = temp.data
            pd = p.data
            sd[salt.mask] = np.nan
            td[salt.mask] = np.nan
            pd[salt.mask] = np.nan
            rho = seawater.dens(sd, td, pd)
            rho = np.ma.masked_where(salt.mask, rho)
            DZ = np.diff(z_w, axis=0)
            bp_arr[tt,:,:] = (g * rho * DZ).sum(axis=0)
        ds1.close()
    bp_mean = np.mean(bp_arr, axis=0)
    bp_anom = bp_arr - bp_mean

# initialize output Dataset
ds2 = nc.Dataset(out_fn, 'w')

# lists of variables to process
if model_type == 'LiveOcean':
    dlist = ['xi_rho', 'eta_rho', 'xi_psi', 'eta_psi', 'ocean_time']
    vn_list2 = [ 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi', 'mask_rho', 'h']
elif model_type == 'Kurapov':
    dlist = ['xi_rho', 'eta_rho', 'ocean_time']
    vn_list2 = [ 'lon_rho', 'lat_rho', 'mask_rho', 'h']
vn_list2t = ['zeta', 'ocean_time']
vn_list3t = ['salt', 'temp']

ds1 = nc.Dataset(fn_list[0])

# Copy dimensions
for dname, the_dim in ds1.dimensions.items():
    if dname in dlist:
        ds2.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

# Copy variables
if model_type == 'Kurapov':
    dsg = nc.Dataset(fng)
for vn in vn_list2:
    if model_type == 'LiveOcean':
        varin = ds1[vn]
    elif model_type == 'Kurapov':
        varin = dsg[vn]
    vv = ds2.createVariable(vn, varin.dtype, varin.dimensions)
    vv.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    if model_type == 'LiveOcean':
        vv[:] = ds1[vn][:]
    elif model_type == 'Kurapov':
        vv[:] = dsg[vn][1030-1:1521, 375-1:615]
if model_type == 'Kurapov':
    dsg.close()
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

if False:
    zfun.ncd(out_fn)
    
if False:
    import matplotlib.pyplot as plt
    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ds = nc.Dataset(out_fn)
    vn = 'bpa'
    pch = ax.pcolormesh(ds['lon_rho'][:], ds['lat_rho'][:], ds[vn][-1,:,:].squeeze())
    fig.colorbar(pch)
    ax.set_title(ds[vn].long_name + ' (' + ds[vn].units + ')')
    pfun.dar(ax)
    plt.show()


