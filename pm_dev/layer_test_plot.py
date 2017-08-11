#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 16:03:28 2017

@author: PM5

Code to test, graphically, the results of layer_extractor.py

"""

# setup
import os
import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
from datetime import datetime, timedelta
import netCDF4 as nc

plp = os.path.abspath('../../LiveOcean/plotting')
if plp not in sys.path:
    sys.path.append(plp)
import pfun

import matplotlib.pyplot as plt

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
outdir0 = Ldir['parent'] + 'ptools_output/layer/'
outdir = outdir0 + Ldir['gtagex'] + '/'
out_fn = outdir + out_name

ds = nc.Dataset(out_fn)

ot = ds['ocean_time'][:]
otu = ds['ocean_time'].units

vn_list2 = [ 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi', 'mask_rho', 'h']
vn_list2t = ['zeta']
vn_list3t = ['salt', 'temp']
#vn_list2t = ['Uwind', 'Vwind', 'zeta']
#vn_list3t = ['salt', 'temp', 'NO3', 'phytoplankton',
#           'zooplankton', 'oxygen', 'TIC', 'alkalinity', 'PH', 'ARAG']

G = dict()
for vn in vn_list2:
    G[vn] = ds[vn][:]
    
# plotting

plt.close('all')

full_list = vn_list2t + vn_list3t
if do_vel:
    full_list = full_list  + ['u', 'v']
    
for vn in full_list:
    fig = plt.figure(figsize=(12,8))
    nplot = 1
    for tlev in [0, -1]:
        ax = fig.add_subplot(1,2,nplot)
        cs = ax.pcolormesh(G['lon_psi'], G['lat_psi'], ds[vn][tlev, 1:-1, 1:-1],
                           cmap='rainbow')
        try:
            tun = ds[vn].units
        except AttributeError:
            tun = ''
        ax.set_title(vn + ' (' + tun + ')')
        fig.colorbar(cs)
        pfun.dar(ax)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        t = ot[tlev]
        fs = 12
        dt = datetime(1970,1,1,0,0,0) + timedelta(days=t/86400)
        ax.text(.95, .075, dt.strftime('%Y-%m-%d'),
            horizontalalignment='right' , verticalalignment='bottom',
            transform=ax.transAxes, fontsize=fs)
        ax.text(.95, .065, dt.strftime('%H:%M') + ' UTC',
            horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes, fontsize=fs)
        nplot += 1
        
plt.show()
