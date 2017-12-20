# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:03:12 2016

@author: PM5
"""

"""
Plot results of tracker using .nc files.

"""

# setup
import os
import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import matfun
import pickle
import numpy as np
import pandas as pd
from datetime import datetime, timedelta

import netCDF4 as nc

plp = os.path.abspath('../../LiveOcean/plotting')
if plp not in sys.path:
    sys.path.append(plp)
import pfun

Ldir = Lfun.Lstart()

if Ldir['env'] == 'pm_mac': # mac version
    pass
elif Ldir['env'] == 'fjord': # fjord version
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt

indir = odir00 = Ldir['parent'] + 'ptools_output/rockfish/'
datadir = '/data1/bbartos/LiveOcean_data/tracker/'

# create the list of run files
m_list_raw = os.listdir(indir)
m_list = []
for m in m_list_raw:
    if ('.nc' in m) and ('grid' not in m):
        m_list.append(m)
m_list.sort()
Npt = len(m_list)
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
my_ndt = int(input('-- Input number (99 for all) -- '))
if my_ndt == 99:
    pass
else:
    m_list = [m_list[my_ndt],]

# retrieve experimental data
exdf = pd.read_csv(datadir + 'rockfish_latlon.csv', index_col = 0)

for inname in m_list:
    
    # compile list of day files
    ds = nc.Dataset(indir + inname)
    dsg = nc.Dataset(indir + 'grid.nc')
    
    # PLOTTING
    
    fig = plt.figure(figsize=(16,8))
    plt.suptitle(species+' rockfish released from '+location+' (2006)')
    ax = fig.add_subplot(1,2,1)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    
    pfun.add_coast(ax)
    aa = [-125, -122, 47, 49]
    ax.axis(aa)
    pfun.dar(ax)
    ax.grid()
    
    # tracks
    ax.plot(ds['lon'][:], ds['lat'][:], '-r', linewidth=0.5, alpha=0.5)
    # ending points
    ax.plot(ds['lon'][-1,:],ds['lat'][-1,:],'ob', markersize=4, label='End')
    # starting points
    ax.plot(ds['lon'][0,:], ds['lat'][0,:], '^y', markersize=10, label='Start')
    ax.legend()
        
    # TIME SERIES
    tdays = (ds['ot'] - ds['ot'][0])/86400.
    # Depth
    ax = fig.add_subplot(3,2,2)
    ax.plot(tdays, ds['z'],'-', alpha=0.25)
    ax.set_ylabel('Z (m)')

    # Distance from Start
    dis = np.zeros(ds['z'].shape)
    for tind in np.arange(1,len(tdays)):
    # change in lat/lon and total distance=sqrt(lat^2+lon^2)
        lat_dis = (ds['lat'][tind,:] - ds['lat'][tind-1,:]) * 111
        lon_dis = ((ds['lon'][tind,:] - ds['lon'][tind-1,:]) * 111 *
                        np.cos(np.nanmean(ds['lat'][tind,:])*np.pi/180))
        dis[tind,:] = dis[tind-1,:] + np.sqrt(lat_dis**2 + lon_dis**2)
    # remove dead larvae
    dis_alive = np.where(np.isfinite(ds['z']), dis, np.nan)
    ax = fig.add_subplot(3,2,4)
    ax.plot(tdays, dis_alive, '-', alpha=0.25)
    ax.set_ylabel('Distance Traveled (km)')
    ax.grid()
    
    # Net Distance from Start
    ndis = np.zeros(ds['z'].shape)
    for tind in np.arange(1,len(tdays)):
    # change in lat/lon and total distance=sqrt(lat^2+lon^2)
        lat_ndis = (ds['lat'][tind,:] - ds['lat'][0,:]) * 111
        lon_ndis = ((ds['lon'][tind,:] - ds['lon'][0,:]) * 111 *
                        np.cos(np.nanmean(ds['lat'][tind,:])*np.pi/180))
        ndis[tind,:] = np.sqrt(lat_ndis**2 + lon_ndis**2)
    # remove dead larvae
    ndis_alive = np.where(np.isfinite(ds['z']), ndis, np.nan)
    ax = fig.add_subplot(3,2,6)
    ax.plot(tdays, ndis_alive, '-', alpha=0.25)
    ax.set_ylabel('Distance From Release(km)')
    ax.set_xlabel('Days')
    ax.grid()

    # save figures
    outfn = indir + inname + '.png'
    plt.savefig(outfn)
    plt.close('all')

