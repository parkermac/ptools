#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 14:35:03 2017

@author: PM5

Code the start looking at the files for the Wendy Schmidt project.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime

import os
import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zfun
import zrfun

alp = os.path.abspath('../../LiveOcean/plotting')
if alp not in sys.path:
    sys.path.append(alp)
import pfun

Ldir = Lfun.Lstart()
indir = Ldir['data'] + 'wendy_schmidt_2007/'

fn = indir + 'usw42_avg_Y2007D181-185.nc'
fng = indir + 'usw42_grd_Z_DZ.nc'

ds = nc.Dataset(fn)
dsg = nc.Dataset(fng)

a = [v for v in ds.variables]
b = [v for v in dsg.variables]

#%% plotting
plt.close('all')

lon = dsg['lon_psi'][:]
lat = dsg['lat_psi'][:]

fig = plt.figure(figsize=(13,8))

vn_list = ['DIATCHL','DIC','NO3','salt','temp','O2']

nr = 2
nc = 3
counter = 0
for vn in vn_list:
    counter += 1
    ax = fig.add_subplot(nr, nc, counter)    
    c = ds[vn][0, -1, :, :].squeeze()    
    cs = ax.pcolormesh(lon, lat, c[1:-1, 1:-1],
                       cmap='rainbow')
    fig.colorbar(cs)
    ax.set_title(vn)
    pfun.dar(ax)

plt.show()
