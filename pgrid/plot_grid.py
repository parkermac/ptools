# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:07:01 2016

@author: PM5
"""

from importlib import reload
import gfun; reload(gfun)
G = gfun.gstart()

import matplotlib.pyplot as plt
import netCDF4 as nc

#import os
import zfun

import numpy as np

#%% select grid file
fn = gfun.select_file(G)
in_fn = G['gdir'] + fn

#%% load the data

ds = nc.Dataset(in_fn)
z = -ds.variables['h'][:]
mask_rho = ds.variables['mask_rho'][:]
plon = ds.variables['lon_psi_ex'][:]
plat = ds.variables['lat_psi_ex'][:]

show_grids = False
NC = 1
if show_grids:
    NC = 2
    lon_dict = dict()
    lat_dict = dict()
    tag_list = ['rho', 'u', 'v', 'psi']
    for tag in tag_list:
        lon_dict[tag] = ds.variables['lon_'+tag][:]
        lat_dict[tag] = ds.variables['lat_'+tag][:]

ds.close()

#%% plotting

cmat = gfun.get_coast()
ax_lims = (plon[0,0], plon[0,-1], plat[0,0], plat[-1,0])
zm = np.ma.masked_where(mask_rho==0, z)

plt.close()

fig = plt.figure(figsize=(10*NC,10))

ax = fig.add_subplot(1,NC,1)
cmap1 = plt.get_cmap(name='terrain')
cs = ax.pcolormesh(plon, plat, zm,
                   vmin=-200, vmax=200, cmap = cmap1)
ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
zfun.dar(ax)
fig.colorbar(cs, ax=ax, extend='both')
ax.set_xlim(ax_lims[:2])
ax.set_ylim(ax_lims[-2:])
ax.set_title(G['gridname'] + '/' + fn)

if show_grids:
    marker_dict = {'rho': 'ok',
                 'u': '>r',
                 'v': '^b',
                 'psi': 'xg'}
    ax = fig.add_subplot(1,NC,2)
    for tag in tag_list:
        if tag == 'rho':
            ax.plot(lon_dict[tag][mask_rho==1], lat_dict[tag][mask_rho==1],
                    marker_dict[tag])
        else:
            ax.plot(lon_dict[tag], lat_dict[tag], marker_dict[tag])
    ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
    zfun.dar(ax)
    ax.set_xlim(ax_lims[:2])
    ax.set_ylim(ax_lims[-2:])

plt.show()
