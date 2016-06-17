# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 16:43:32 2016

@author: PM5

This plots river tacks.
"""

dir0 = '/Users/PM5/Documents/'

import os
import sys
alp = os.path.abspath(dir0 + 'LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()
import zfun
import matfun

import netCDF4 as nc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%% get river info
gname = 'test'
outdir0 = Ldir['LOo'] + 'grids/'
outdir = outdir0 + gname +'/'
trackdir = outdir + 'tracks/'

fn_ri = outdir + 'river_info.csv'

df = pd.read_csv(fn_ri, index_col='rname')

#%% get a bathymetry file

dir0 = '/Users/PM5/Documents/'
outdir = dir0 +'ptools_output/pgrid/'
gname_base = 'test'
gname_tag = '_m02_r00_s00_x00'
gname = gname_base + gname_tag + '.nc'
g_fn = outdir + gname

ds = nc.Dataset(g_fn)
z = ds.variables['z'][:]
mask_rho = ds.variables['mask_rho'][:]
plon = ds.variables['lon_psi_ex'][:]
plat = ds.variables['lat_psi_ex'][:]
ds.close()

ax_lims = (plon[0,0], plon[0,-1], plat[0,0], plat[-1,0])

#%% try carving river channels

zm = mask_rho==0

plon_vec = plon[0,:].flatten()
plat_vec = plat[:,0].flatten()

for rn in df.index:
    depth = df.ix[rn, 'depth']
    fn_tr = trackdir + rn + '_track.csv'
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
        zm[jj0, ii0] = False
        # and this sets the depth in those places, if needed,
        # "carving the river channel"
        for ff in range(len(ii0)):
            JJ = jj0[ff]
            II = ii0[ff]
            z[JJ, II] = np.min((z[JJ, II], -depth))

# make a masked array for plotting
zma = np.ma.masked_where(zm, z)

#%% coast

c_dir = dir0 + 'tools_data/geo_data/coast/'
c_file = 'pnw_coast_combined.mat'
c_fn = c_dir + c_file
cmat = matfun.loadmat(c_fn)

#%% plotting

plt.close()

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
cmap1 = plt.get_cmap(name='rainbow')
cs = ax.pcolormesh(plon, plat,zma,
                   vmin=-200, vmax=200, cmap = cmap1)
ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
zfun.dar(ax)
fig.colorbar(cs, ax=ax, extend='both')
ax.set_xlim(ax_lims[:2])
ax.set_ylim(ax_lims[-2:])
ax.set_title(gname)

for rn in df.index:
    fn_tr = trackdir + rn + '_track.csv'
    df_tr = pd.read_csv(fn_tr, index_col='ind')
    x = df_tr['lon'].values
    y = df_tr['lat'].values
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[-1], y[-1], '*r')

plt.show()