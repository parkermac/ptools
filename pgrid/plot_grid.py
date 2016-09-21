# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:07:01 2016

@author: PM5
"""

from importlib import reload
import gfun; reload(gfun)
G = gfun.gstart()
import pfun

import matplotlib.pyplot as plt
import netCDF4 as nc

import numpy as np

#%% select grid file
LO=False
if LO==True:
    fn = gfun.select_file(G, LO=True)
    in_fn = fn
elif LO==False:
    fn = gfun.select_file(G)
    in_fn = G['gdir'] + fn

#%% load the data

ds = nc.Dataset(in_fn)
z = -ds.variables['h'][:]
mask_rho = ds.variables['mask_rho'][:]
if LO==True:
    # because older grids did not have lon,lat_psi_ex we create this
    # as an extension of lon,lat_psi
    plon0 = ds.variables['lon_psi'][:]
    plat0 = ds.variables['lat_psi'][:]
    dx = plon0[0,1] - plon0[0,0]
    dy = plat0[1,0] - plat0[0,0]
    ny0, nx0 = plon0.shape
    plon = np.nan * np.ones((ny0+2, nx0+2))
    plat = np.nan * np.ones((ny0+2, nx0+2))
    plon[1:-1, 1:-1] = plon0
    plat[1:-1, 1:-1] = plat0
    plon[:,0] = plon0[0,0] - dx
    plon[:,-1] = plon0[0,-1] + dx
    plon[0,:] = plon[1,:]
    plon[-1,:] = plon[-2,:]
    plat[0,:] = plat0[0,0] - dy
    plat[-1,:] = plat0[-1,0] + dy
    plat[:,0] = plat[:,1]
    plat[:,-1] = plat[:,-2]

elif LO==False:
    plon = ds.variables['lon_psi_ex'][:]
    plat = ds.variables['lat_psi_ex'][:]

show_grids = False
NC = 1
if show_grids:
    NC = 2
    lon_dict = dict()
    lat_dict = dict()
    mask_dict = dict()
    tag_list = ['rho', 'u', 'v', 'psi']
    for tag in tag_list:
        lon_dict[tag] = ds.variables['lon_'+tag][:]
        lat_dict[tag] = ds.variables['lat_'+tag][:]
        mask_dict[tag] = ds.variables['mask_'+tag][:]


#%% plotting

ax_lims = (plon[0,0], plon[0,-1], plat[0,0], plat[-1,0])
zm = np.ma.masked_where(mask_rho==0, z)

plt.close()

fig = plt.figure(figsize=(10*NC,10))

ax = fig.add_subplot(1,NC,1)
cmap1 = plt.get_cmap(name='terrain')
cs = ax.pcolormesh(plon, plat, zm,
                   vmin=-200, vmax=200, cmap = cmap1)
fig.colorbar(cs, ax=ax, extend='both')
pfun.add_coast(ax)
pfun.dar(ax)
ax.axis(pfun.get_aa(ds))
ax.set_title(G['gridname'] + '/' + fn)

(rowmax, colmax) = np.unravel_index(np.argmax(zm), zm.shape)
zmax = zm[rowmax, colmax]
print('Max z = ' + str(zmax))
lon_rho = ds['lon_rho'][:]
lat_rho = ds['lat_rho'][:]
ax.plot(lon_rho[rowmax, colmax], lat_rho[rowmax, colmax], '*m', markersize=20)

if show_grids:
    marker_dict = {'rho': 'ok',
                 'u': '>r',
                 'v': '^b',
                 'psi': 'xg'}
    ax = fig.add_subplot(1,NC,2)
    for tag in tag_list:
        ax.plot(lon_dict[tag][mask_dict[tag]==1], lat_dict[tag][mask_dict[tag]==1],
                marker_dict[tag])
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(pfun.get_aa(ds))

ds.close()

plt.show()
