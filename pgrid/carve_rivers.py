# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 16:43:32 2016

@author: PM5

This plots river tacks.
"""

from importlib import reload
import gfun; reload(gfun)
G = gfun.gstart()

import zfun

import os
import shutil
import netCDF4 as nc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seawater as sw

#%% get river info
ri_name = 'sng_2016_06'

ri_dir0 = G['dir0'] + 'ptools_output/river/'
ri_dir = ri_dir0 + ri_name +'/'
rit_dir = ri_dir + 'tracks/'

ri_fn = ri_dir + 'river_info.csv'

df = pd.read_csv(ri_fn, index_col='rname')

#%% select grid file
fn = gfun.select_file(G)
in_fn = G['gdir'] + fn
# create new file name
fn_new = gfun.increment_filename(fn, tag='_r')
out_fn = G['gdir'] + fn_new

#%% get the grid data

ds = nc.Dataset(in_fn)
z = ds.variables['z'][:]
mask_rho = ds.variables['mask_rho'][:]
plon = ds.variables['lon_psi_ex'][:]
plat = ds.variables['lat_psi_ex'][:]
ds.close()

ax_lims = (plon[0,0], plon[0,-1], plat[0,0], plat[-1,0])

#%% try carving river channels

m = mask_rho==0

plon_vec = plon[0,:].flatten()
plat_vec = plat[:,0].flatten()

ilon_dict = dict()
ilat_dict = dict()
dir_dict = dict()

for rn in df.index:
    depth = df.ix[rn, 'depth']
    fn_tr = rit_dir + rn + '.csv'
    df_tr = pd.read_csv(fn_tr, index_col='ind')
    x = df_tr['lon'].values
    y = df_tr['lat'].values
    # This unmasks it in the places where the
    # river crosses a tile
    for I in range(len(x)-1):
        xx = np.linspace(x[I], x[I+1], 10)
        yy = np.linspace(y[I], y[I+1], 10)
        ii0, ii1, ifr = zfun.get_interpolant(xx, plon_vec, extrap_nan=True)
        jj0, jj1, jfr = zfun.get_interpolant(yy, plat_vec, extrap_nan=True)
        # drop extrapolated points
        ii0 = ii0[~np.isnan(ifr) & ~np.isnan(jfr)]
        jj0 = jj0[~np.isnan(ifr) & ~np.isnan(jfr)]
        # this unmasks points crossed by the river
        m[jj0, ii0] = False
        # and this sets the depth in those places, if needed,
        # "carving the river channel"
        for ff in range(len(ii0)):
            JJ = jj0[ff]
            II = ii0[ff]
            z[JJ, II] = np.min((z[JJ, II], -depth))
    # this creates information about the index and direction of the river
    ilon_dict[rn] = II
    ilat_dict[rn] = JJ
    # phaseangle is degrees -180:180 with 0 = East
    dist, phaseangle = sw.dist([y[I], y[I+1]], [x[I], x[I+1]])
    dir_dict[rn] = phaseangle

#%% figure out river locations
ji_dict = dict()
for rn in df.index:
    ii = ilon_dict[rn]
    jj = ilat_dict[rn]
    ji = np.array([jj,ii])
    ph = dir_dict[rn]
    if -45 <= ph and ph <= 45:
        JI = np.array([0,1])
    elif 135 <= ph and ph < 45:
        JI = np.array([1,0])
    elif -135 <= ph and ph < -45:
        JI = np.array([-1,0])
    elif np.abs(ph) > 135:
        JI = np.array([0,-1])
    is_right = False
    while is_right == False:
        ji = ji + JI
        if mask_rho[ji[0],ji[1]] == 0:
            ji_dict[rn] = ji
            print(rn)
            is_right = True

#%% plotting

zm = np.ma.masked_where(m, z)

cmat = gfun.get_coast()

plt.close()

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
cmap1 = plt.get_cmap(name='rainbow')
cs = ax.pcolormesh(plon, plat,zm,
                   vmin=-200, vmax=200, cmap = cmap1)
ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
zfun.dar(ax)
fig.colorbar(cs, ax=ax, extend='both')
ax.set_xlim(ax_lims[:2])
ax.set_ylim(ax_lims[-2:])
ax.set_title(G['gridname'])

for rn in df.index:
    fn_tr = rit_dir + rn + '.csv'
    df_tr = pd.read_csv(fn_tr, index_col='ind')
    x = df_tr['lon'].values
    y = df_tr['lat'].values
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[-1], y[-1], '*r')

plt.show()

#%% Save the output file

# create the new mask
mask_rho[m == True] = 0
mask_rho[m == False] = 1

print('Creating ' + out_fn)
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist
shutil.copyfile(in_fn, out_fn)
ds = nc.Dataset(out_fn, 'a')
ds['mask_rho'][:] = mask_rho
ds['z'][:] = z
ds.close()
