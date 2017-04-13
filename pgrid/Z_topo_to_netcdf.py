#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 16:18:52 2017

@author: PM5

Code to convert bathymetry files from matlab to netcdf.
"""

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

dir0 = '/Users/PM5/Documents/'

import os
import sys
pth = os.path.abspath(dir0 +'LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import zfun
    
pth = os.path.abspath(dir0 +'LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

from importlib import reload
import gfun_utility as gfu
reload(gfu)

tdir0 = dir0 + 'tools_data/geo_data/topo/'

tdict = {'srtm15': 'topo15.grd',
         'cascadia': 'cascadia_gridded.mat',
         'psdem': 'PS_183m.mat',
         'ttp_patch': 'TTP_Regional_27m_patch.mat'}

#tdict = {'psdem': 'PS_183m.mat'}

plt.close('all')

for k in tdict.keys():
    t_fn = tdir0 + k + '/' + tdict[k]
    print('Opening %s/%s' % (k, tdict[k]))
    if k == 'srtm15':        
        lon_vec = np.array([-135, -122])
        lat_vec = np.array([35, 53])
        lon, lat, z = gfu.load_bathy2(t_fn, lon_vec, lat_vec)
    else:
        lon, lat, z = gfu.load_bathy(t_fn)
        
    z = np.ma.masked_where(np.isnan(z), z)
    vv = int(np.nanstd(z))/2
    aa = [lon[0], lon[-1], lat[0], lat[-1]]
    
    if False:
        # plot from original data
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        plon = lon[:-1] + np.diff(lon)/2
        plat = lat[:-1] + np.diff(lat)/2
        cs = ax.pcolormesh(plon, plat, z[1:-1, 1:-1],
                      cmap='terrain', vmin=-vv, vmax=vv)
        fig.colorbar(cs, ax=ax, extend='both')
        pfun.dar(ax)
        pfun.add_coast(ax)
        ax.axis(aa)
        R = zfun.earth_rad(np.mean(lat))
        dx = R * np.cos(np.pi*lat/180) * (np.pi*np.diff(lon[:2])/180)
        dy = R * (np.pi*np.diff(lat)/180)
        minres = np.min([np.min(dx), np.min(dy)])
        ax.set_title('%s/%s: minres = %d m' % (k, tdict[k], int(minres)))
        plt.show()
    
    # write the data to NetCDF
    out_fn = tdir0 + k + '/' + tdict[k].strip('.mat') + '.nc'
    foo = nc.Dataset(out_fn, 'w')
    L = len(lon)
    M = len(lat)
    foo.createDimension('x', L)
    foo.createDimension('y', M)
    lon_var = foo.createVariable('lon', float, ('x'))
    lat_var = foo.createVariable('lat', float, ('y'))
    z_var = foo.createVariable('z', float, ('y', 'x'))
    lon_var[:] = lon
    lat_var[:] = lat
    z_var[:] = z
    foo.close()

    if True:
        # plot from NetCDF
        ds = nc.Dataset(out_fn)
        lon = ds['lon'][:]
        lat = ds['lat'][:]
        z = ds['z'][:]
        ds.close()
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        plon = lon[:-1] + np.diff(lon)/2
        plat = lat[:-1] + np.diff(lat)/2
        cs = ax.pcolormesh(plon, plat, z[1:-1, 1:-1],
                      cmap='terrain', vmin=-vv, vmax=vv)
        fig.colorbar(cs, ax=ax, extend='both')
        pfun.dar(ax)
        pfun.add_coast(ax)
        ax.axis(aa)
        R = zfun.earth_rad(np.mean(lat))
        dx = R * np.cos(np.pi*lat/180) * (np.pi*np.diff(lon[:2])/180)
        dy = R * (np.pi*np.diff(lat)/180)
        minres = np.min([np.min(dx), np.min(dy)])
        ax.set_title('%s/%s: minres = %d m' % (k, tdict[k], int(minres)))
        plt.show()
