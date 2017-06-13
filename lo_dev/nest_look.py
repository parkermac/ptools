#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 30 09:23:03 2017

@author: PM5

Code to explore some ROMS output.  Specifically for debugging
my new nesting run.
"""

#%% setup
import os
import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
from importlib import reload
import Lfun
Ldir = Lfun.Lstart(gridname='sal0', tag='f1')
Ldir['gtagex'] = Ldir['gtag'] + '_' + 'r820_nest'
import zfun
import zrfun

alp = os.path.abspath('../../LiveOcean/plotting')
if alp not in sys.path:
    sys.path.append(alp)
import pfun

#import roms_plots; reload(roms_plots)

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

#%%
in_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f2013.01.01/'
in_fn = 'ocean_his_0001.nc'
fn = in_dir + in_fn
G, S, T = zrfun.get_basic_info(fn)
ds = nc.Dataset(fn)
zeta = ds['zeta'][:].squeeze()
ds.close()

#%%
in_fn = 'ocean_rst.nc'
fn = in_dir + in_fn
ds = nc.Dataset(fn)
print(50*'*')
a = [print(ds[v]) for v in ds.variables]
print('\n\n')

t = ds['ocean_time'][:]
u = ds['u'][0,1,-1,:,:].squeeze()
v = ds['v'][0,1,-1,:,:].squeeze()
#ds.close()


#%% derived quantities
h = G['h']
hm = np.ma.masked_where(zeta.mask==True, h)

th = zeta+hm

#%% PLOTTING

plt.close('all')

fig = plt.figure(figsize=(13,8))
ax = fig.add_subplot(111)
#cs = ax.pcolormesh(G['lon_psi'], G['lat_psi'], th[1:-1, 1:-1])
cs = ax.pcolormesh(G['lon_u'], G['lat_u'], u)
pfun.dar(ax)
fig.colorbar(cs)

plt.show()

