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
z_alt = -ds.variables['h_alt'][:]
mask_rho = ds.variables['mask_rho'][:]
mask_rho_alt = ds.variables['mask_rho_alt'][:]

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

NC = 1 # first guess at number of columns for plot

show_grids = False
if show_grids:
    NC += 1
    lon_dict = dict()
    lat_dict = dict()
    mask_dict = dict()
    tag_list = ['rho', 'u', 'v', 'psi']
    for tag in tag_list:
        lon_dict[tag] = ds.variables['lon_'+tag][:]
        lat_dict[tag] = ds.variables['lat_'+tag][:]
        mask_dict[tag] = ds.variables['mask_'+tag][:]
        
show_sections = True
if show_sections:
    NC += 1

#%% plotting

ax_lims = (plon[0,0], plon[0,-1], plat[0,0], plat[-1,0])

zm = np.ma.masked_where(mask_rho_alt == 0, z_alt)

plt.close()

fig = plt.figure(figsize=(10*NC,10))

ax1 = fig.add_subplot(1,NC,1)
cmap1 = plt.get_cmap(name='rainbow') # terrain, viridis
cs = ax1.pcolormesh(plon, plat, zm,
                   vmin=-500, vmax=500, cmap = cmap1)
fig.colorbar(cs, ax=ax1, extend='both')
pfun.add_coast(ax1)
pfun.dar(ax1)
ax1.axis(pfun.get_aa(ds))
ax1.set_title(G['gridname'] + '/' + fn)

if False:
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
    
if show_sections:
    color_list = ['orange', 'gold', 'greenyellow', 'lightgreen',
                    'aquamarine', 'cadetblue', 'royalblue', 'purple']
    lon_rho = ds['lon_rho'][:]
    lat_rho = ds['lat_rho'][:]
    NS = 5 # number of sections
    for ss in range(NS):
        ax = fig.add_subplot(NS,NC,2*NS - 2*(ss+1) + 2)
        x = lon_rho[0, :]
        jj = int(lon_rho.shape[0]/(NS+1) * (ss+1))
        y = z[jj, :]/100
        y_alt = z_alt[jj, :]/100
        ax.plot(x, y, '-k', x, y_alt, '-r')
        ax.plot(x, 0*x, '-', color=color_list[ss], linewidth=2)
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(-3, 1)
        if ss == NS-1:
            ax.set_title('Z/100 (red = alt)')
        
        ax1.plot([x[0], x[-1]], [lat_rho[jj,0], lat_rho[jj, -1]], '-',
            color=color_list[ss], linewidth=2)

ds.close()

plt.show()
pfun.topfig()
