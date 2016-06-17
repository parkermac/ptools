# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:07:01 2016

@author: PM5
"""

from importlib import reload
import gfun; reload(gfun)
gridname, dir0, pgdir, gdir = gfun.gstart()

import matplotlib.pyplot as plt
import netCDF4 as nc

import os
import zfun

import numpy as np

#%% load grid file

if True:
    # interactive selection
    print('\n%s\n' % '** Choose file to edit **')
    fn_list_raw = os.listdir(gdir)
    fn_list = []
    for item in fn_list_raw:
        if item[-3:] == '.nc':
            fn_list.append(item)
    Nfn = len(fn_list)
    fn_dict = dict(zip(range(Nfn), fn_list))
    for nfn in range(Nfn):
        print(str(nfn) + ': ' + fn_list[nfn])
    my_nfn = int(input('-- Input number -- '))
    fn = fn_dict[my_nfn]
else:
    fn = 'grid_m00_r00_s00_x00.nc'

in_fn = gdir + fn

#%% coastline

cmat = gfun.get_coast()

#%% Test: retrieve the output

ds = nc.Dataset(in_fn)
z = ds.variables['z'][:]
m = ds.variables['mask_rho'][:]
zm = np.ma.masked_where(m==0, z)
plon = ds.variables['lon_psi_ex'][:]
plat = ds.variables['lat_psi_ex'][:]
ds.close()

ax_lims = (plon[0,0], plon[0,-1], plat[0,0], plat[-1,0])

#%% plotting

plt.close()

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
cmap1 = plt.get_cmap(name='terrain')
cs = ax.pcolormesh(plon, plat, zm,
                   vmin=-200, vmax=200, cmap = cmap1)
ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
zfun.dar(ax)
fig.colorbar(cs, ax=ax, extend='both')
ax.set_xlim(ax_lims[:2])
ax.set_ylim(ax_lims[-2:])
ax.set_title(gridname + '/' + fn)

plt.show()