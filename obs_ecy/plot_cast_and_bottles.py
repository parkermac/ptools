#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots data from the Department of Ecology, combining data from CTD casts
and bottles at the same station.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()
import zfun

# +++ load ecology CTD cast and bottle data +++
dir0 = Ldir['parent'] + 'ptools_data/ecology/'
# load processed station info and data
sta_df = pd.read_pickle(dir0 + 'sta_df.p')
year = 2017
Casts = pd.read_pickle(dir0 + 'Casts_' + str(year) + '.p')
Bottles = pd.read_pickle(dir0 + 'Bottles_' + str(year) + '.p')

# sta_to_plot = ['PSB003']
sta_to_plot = [sn for sn in sta_df.index if (('GYS' not in sn) and ('WPA' not in sn))]

Bc = pd.DataFrame()
for station in sta_to_plot:
    
    # loop over all casts at this station
    print(' - Working on: ' + station)
    
    casts = Casts[Casts['Station'] == station]   
    casts = casts.set_index('Date')    
    # identify a single cast by its DATE
    calldates = casts.index
    castdates = calldates.unique() # a short list of unique dates (1 per cast)
    
    bottles = Bottles[Bottles['Station'] == station]   
    bottles = bottles.set_index('Date')    
    # identify a single cast by its DATE
    balldates = bottles.index
    bottledates = balldates.unique() # a short list of unique dates (1 per cast)
    
    dates = bottledates.intersection(castdates)
    
    for dd in dates:
        # NOTE the brackets around [dd] keep the result as a DataFrame even if
        # we are only pulling out a single row.
        ca = casts.loc[[dd],:]
        bo = bottles.loc[[dd],:]
        ca = ca.set_index('Z')
        bo = bo.set_index('Z')
        bc = bo.copy() # a DataFrame to store both bottles and casts, but only at shared depths
        
        vn_list = ['Salinity', 'Temperature', 'Sigma', 'Chl', 'DO', 'Turb']
        for vn in vn_list:
            bc[vn] = np.nan
        
        for Z in bo.index:
            i0, i1, fr = zfun.get_interpolant(np.array([Z]), ca.index.values, extrap_nan=False)
            for vn in vn_list:
                cv = ca[vn]
                bc.loc[Z,vn] = (1-fr)*cv.iloc[int(i0)] + fr*cv.iloc[int(i1)]
                if False: # test of the interpolation
                    exp = Z
                    obs = (1-fr)*cv.index[int(i0)] + fr*cv.index[int(i1)]
                    print('obs = %0.2f exp = %0.2f' % (obs, exp))
                    # result: looks good - only different at ends
        bc['Date'] = dd
        bc['Z'] = bc.index
        Bc = pd.concat((Bc,bc),ignore_index=True, sort=False)

# plotting
plt.close('all')
fig = plt.figure(figsize=(13,8))

# seasonal values
bcy = Bc.set_index('Date')
bcy = bcy[bcy['Znom']==-30]
winter_mask = (bcy.index.month==11) | (bcy.index.month==12) | (bcy.index.month==1) | (bcy.index.month==2)
summer_mask = (bcy.index.month==5) | (bcy.index.month==6) | (bcy.index.month==7) | (bcy.index.month==8)
bcw = bcy[winter_mask]
bcs = bcy[summer_mask]

vn_list = ['Temperature', 'Chl', 'DIN', 'DO']

nplt = 1
NR = 2
NC = 2
for vn in vn_list:
    ax = fig.add_subplot(NR, NC, nplt)
    ax.plot(bcs['Salinity'].values, bcs[vn].values, '*r', label='Summer', alpha=0.5)
    ax.plot(bcw['Salinity'].values, bcw[vn].values, 'ob', label='Winter', alpha=0.5)
    ax.set_xlabel('Salinity')
    ax.set_ylabel(vn)
    ax.grid(True)
    if nplt==1:
        ax.legend()
    nplt += 1
fig.suptitle('Puget Sound ' + str(year) + ' Ecology Data at 30 m Depth')


plt.show()
