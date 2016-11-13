# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 14:31:46 2016

@author: PM5

Code to check on the time variability of low-passed surface height.
"""

#%% Imports

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()
import zrfun
import zfun

from datetime import datetime, timedelta
import netCDF4 as nc
import matplotlib.pyplot as plt

#%% get zeta time series

R_in_dir0 = Ldir['parent'] + 'roms/output/salish_2006_4_lp/'

dt0 = datetime(2006,7,1)

lon0 = -123
lat0 = 48.25

counter = 0
eta_list = []
dt_list = []
for ndays in range(31): # 209 is 2006.07.29

    dt = dt0 + timedelta(days=ndays)

    f_string = 'f' + dt.strftime('%Y.%m.%d')
    R_in_dir = R_in_dir0 + f_string + '/'
    R_fn = R_in_dir + 'low_passed.nc'
    ds = nc.Dataset(R_fn)

    if counter == 0:
        [G, S] = zrfun.get_basic_info(R_fn, getT=False)
        lon = G['lon_rho'][0,:]
        lat = G['lat_rho'][:,0]
        ii = zfun.find_nearest_ind(lon,lon0)
        jj = zfun.find_nearest_ind(lat,lat0)

    [T] = zrfun.get_basic_info(R_fn, getG=False, getS=False)
    eta_list.append(ds['zeta'][0, jj, ii])
    dt_list.append(T['tm'])

    ds.close()
    counter += 1

#%% plotting

import pandas as pd

es = pd.Series(dict(zip(dt_list, eta_list)))

es.plot()

plt.show()

