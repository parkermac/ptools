# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:07:01 2016

@author: PM5
"""

from importlib import reload
import gfun; reload(gfun)
G = gfun.gstart()
import pfun; reload(pfun)

import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd
import numpy as np

# select grid file

using_old_grid=False
# Set this to True to look at grids we have already created,
# e.g. ones currently in use for LiveOcean.
# Set it to False when interacting with grids from pgrid_output.

if using_old_grid==True:
    fn = gfun.select_file(G, using_old_grid=True)
    in_fn = fn
elif using_old_grid==False:
    fn = gfun.select_file(G)
    in_fn = G['gdir'] + fn

# load the data
ds = nc.Dataset(in_fn)
z = -ds.variables['h'][:]
mask_rho = ds.variables['mask_rho'][:]

plon, plat = gfun.get_plon_plat(using_old_grid, ds)
ax_lims = (plon[0,0], plon[0,-1], plat[0,0], plat[-1,0])

NC = 1 # first guess at number of columns for plot
   
flag_show_grids = False
if flag_show_grids:
    NC += 1
    ax_grids = fig.add_subplot(1,NC,2)
        
flag_show_sections = True
if flag_show_sections:
    NC += 1

# plotting

zm = np.ma.masked_where(mask_rho == 0, z)

plt.close()

fig = plt.figure(figsize=(10*NC,10))

ax1 = fig.add_subplot(1,NC,1)
cmap1 = plt.get_cmap(name='viridis') # terrain, viridis
cs = ax1.pcolormesh(plon, plat, zm,
                   vmin=-200, vmax=10, cmap = cmap1)
fig.colorbar(cs, ax=ax1, extend='both')
pfun.add_coast(ax1)
pfun.dar(ax1)
ax1.axis(ax_lims)
ax1.set_title(G['gridname'] + '/' + fn)
                    
gfun.add_river_tracks(G, ds, ax1)
   
if flag_show_sections:
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
        ax.plot(x, y, '-r', linewidth=2)
        ax.plot(x, 0*x, '-', color=color_list[ss], linewidth=2)
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(-3, 1)
        if ss == NS-1:
            ax.set_title('Z/(100 m)')
        
        ax1.plot([x[0], x[-1]], [lat_rho[jj,0], lat_rho[jj, -1]], '-',
            color=color_list[ss], linewidth=2)

ds.close()

plt.show()
pfun.topfig()
