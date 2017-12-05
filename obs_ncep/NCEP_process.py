# -*- coding: utf-8 -*-
"""
Code to explore NCEP reanalysis fields.
"""
import netCDF4 as ncf
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pandas as pd
from importlib import reload
import os
import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import zfun
reload(zfun)

indir = '/Users/PM5/Documents/tools_data/obs_data/reanalysis/'

# where are the data csv files
dir0 = '/Users/PM5/Documents/'
indir = dir0 + 'tools_data/obs_data/reanalysis/'

outdir0 = dir0 + 'ptools_output/env_data/'
outdir = outdir0 + 'reanalysis/'

try:
    os.mkdir(outdir0)
except OSError:
    pass # assume OSError was raised because directory already exists
try:
    os.mkdir(outdir)
except OSError:
    pass # assume OSError was raised because directory already exists

# 365 daily maps
#infile = 'dswrf.sfc.day.1981-2010.ltm.nc'

# 814 monthly maps
infile = 'dswrf.sfc.mon.mean.nc'
fn = indir + infile

#zfun.ncd(fn)

# load data
ds = ncf.Dataset(fn)
lon = ds.variables['lon'][:]
lat = ds.variables['lat'][:]
f = ds.variables['dswrf'][:]
f_units = ds.variables['dswrf'].units
t = ds.variables['time'][:]
tu = ds.variables['time'].units
ds.close()

# make a datetime axis
if tu == 'hours since 1800-01-01 00:00:0.0':
    dt = []
    for tt in t:
        dt.append(datetime(1800,1,1) + timedelta(tt/24.))
        
# re-pack latitude from south to north
lat = lat[::-1]
f = f[:, ::-1, :]
       
# get a time series
x = np.array([-122. + 360.,])
y = np.array([47.5,])

xi = zfun.get_interpolant(x, np.array(lon))
yi = zfun.get_interpolant(y, np.array(lat))

ii = int(xi[0, 0])
jj = int(yi[0, 0])
tt = 5

# time series
ff = f[:, jj, ii]

#%% make monthly climatology
fdf = pd.DataFrame(index=dt, columns=['ff', 'anom'])
fdf['ff'] = ff
fdf['anom'] = ff
for mo in range(1,13):
    fdf.ix[fdf.index.month==mo, 'anom'] -= (
        fdf.ix[fdf.index.month==mo, 'ff'].mean(axis=0) )
        
fdf.to_pickle(outdir + 'df_NCEP_dswrf_Seattle.p')
          
#%%  Plotting
fig_size = (20,10)
plt.close()
fig = plt.figure(figsize=fig_size)
ax = fig.add_subplot(121)
cmap = plt.get_cmap(name='jet')
cs = ax.pcolormesh(lon, lat, f[tt,:,:],
    vmin=0, vmax=400,  cmap=cmap)
ax.plot(lon[ii], lat[jj], 'ok', markersize=5)
ax.set_xlim(-150+360, -110+360)
ax.set_ylim(10, 70)
zfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Downward Solar Radiation Flux at surface ' + f_units)
fig.colorbar(cs)

ax = fig.add_subplot(222)
ax.plot(dt, ff)
ax.plot(dt[tt], ff[tt], 'ok', markersize=5)
ax.set_ylim(0, 400)

ax = fig.add_subplot(224)
# fdf['anom'].plot()
fdf['anom'].resample('A', how='mean').plot()
ax.grid()

plt.show()
