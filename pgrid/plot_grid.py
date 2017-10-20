# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:07:01 2016

@author: PM5
"""


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', nargs='?', default='',
        type=str)
args = parser.parse_args()

from importlib import reload
import gfun
reload(gfun)
if len(args.gridname) > 0:
    Gr = gfun.gstart(gridname=args.gridname)
else:
    Gr = gfun.gstart()
import pfun
reload(pfun)
import gfun_plotting as gfp

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

# select grid file

using_old_grid = False
# Set this to True to look at grids we have already created,
# e.g. ones currently in use for LiveOcean.
# Set it to False when interacting with grids from pgrid_output.

if using_old_grid==True:
    fn = gfun.select_file(Gr, using_old_grid=True)
    in_fn = fn
elif using_old_grid==False:
    fn = gfun.select_file(Gr)
    in_fn = Gr['gdir'] + fn
    in_fn0 = Gr['gdir'] + 'grid_m00_r00_s00_x00.nc'
    import pickle
    dch = pickle.load(open(Gr['gdir'] + 'choices.p', 'rb'))

# load the data
ds = nc.Dataset(in_fn)
z = -ds.variables['h'][:]
if using_old_grid==False:
    ds0 = nc.Dataset(in_fn0)
    z0 = -ds0.variables['h'][:]
mask_rho = ds.variables['mask_rho'][:]

plon, plat = gfp.get_plon_plat(using_old_grid, ds)
buff = 0.05*(plat[-1,0]-plat[0,0])
ax_lims = (plon[0,0]-buff, plon[0,-1]+buff, plat[0,0]-buff, plat[-1,0]+buff)

zm = np.ma.masked_where(mask_rho == 0, z)
if using_old_grid==False:
    z0m = np.ma.masked_where(mask_rho == 0, z0)

# plotting
plt.close()

# set number of columns for plot 
NC = 1 # first guess
flag_show_grids = False
if flag_show_grids:
    NC += 1
    icg = NC       
flag_show_sections = True
if flag_show_sections:
    NC += 1
    ics = NC
fig = plt.figure(figsize=(8*NC,8))

#ax_grids = fig.add_subplot(1,NC,2)

ax1 = fig.add_subplot(1,NC,1)
cmap1 = plt.get_cmap(name='terrain') # terrain, viridis
cs = ax1.pcolormesh(plon, plat, zm,
                   vmin=-200, vmax=100, cmap = cmap1)
fig.colorbar(cs, ax=ax1, extend='both')
pfun.add_coast(ax1)
pfun.dar(ax1)
ax1.axis(ax_lims)
ax1.set_title(Gr['gridname'] + '/' + fn)
ax1.text(.95, .05, str(mask_rho.shape), horizontalalignment='right',
         transform=ax1.transAxes)                   
gfp.add_river_tracks(Gr, ds, ax1)
   
if flag_show_sections:
    color_list = ['orange', 'gold', 'greenyellow', 'lightgreen',
                    'aquamarine', 'cadetblue', 'royalblue', 'purple']
    lon_rho = ds['lon_rho'][:]
    lat_rho = ds['lat_rho'][:]
    NS = 5 # number of sections
    for ss in range(NS):
        ax = fig.add_subplot(NS,NC,ics*NS - ics*(ss+1) + ics)
        x = lon_rho[0, :]
        jj = int(lon_rho.shape[0]/(NS+1) * (ss+1))
        y = z[jj, :]/100
        y0 = z0[jj, :]/100
        ax.plot(x, y, '-r', linewidth=1)
        if using_old_grid==False:
            ax.plot(x, y0, '-c', linewidth=1)
        ax.plot(x, 0*x, '-', color=color_list[ss], linewidth=1)
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(-5, 1)
        if ss == NS-1:
            ax.set_title('Z/(100 m)')
        
        ax1.plot([x[0], x[-1]], [lat_rho[jj,0], lat_rho[jj, -1]], '-',
            color=color_list[ss], linewidth=2)
        
if flag_show_grids:
    # NOTE: you need to have run make_extras.py for this to work
    lon_dict = dict()
    lat_dict = dict()
    mask_dict = dict()
    tag_list = ['rho', 'u', 'v', 'psi']
    for tag in tag_list:
        lon_dict[tag] = ds.variables['lon_'+tag][:]
        lat_dict[tag] = ds.variables['lat_'+tag][:]
        mask_dict[tag] = ds.variables['mask_'+tag][:]
    marker_dict = {'rho': 'ok',
                 'u': '>r',
                 'v': '^b',
                 'psi': 'xg'}
    
    ax = fig.add_subplot(1,NC,icg)
    for tag in tag_list:
        ax.plot(lon_dict[tag][mask_dict[tag]==1], lat_dict[tag][mask_dict[tag]==1],
                marker_dict[tag])
        ax.plot(lon_dict[tag][mask_dict[tag]==0], lat_dict[tag][mask_dict[tag]==0],
                marker_dict[tag], alpha = .2)
    pfun.dar(ax)
    ax.set_xlim(ax_lims[:2])
    ax.set_ylim(ax_lims[-2:])
    
    pfun.add_coast(ax)
    ax.axis(ax_lims)
    gfp.add_river_tracks(Gr, ds, ax)

ds.close()
ds0.close()

plt.show()
pfun.topfig()
