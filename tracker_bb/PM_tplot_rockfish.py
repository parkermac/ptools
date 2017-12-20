# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:03:12 2016

@author: PM5
"""

"""
Plot results of tracker.

Files are in places like
/data1/bbartos/LiveOcean_output/tracks/
MoSSea_rockfish_rk4_ndiv1_forward_surfaceFalse_turbTrue_windage0_boundaryreflect/

rockfish_2006_Experiment_3_1/

with names day_00000.p through day_00179.p

Each Experiment folder is about 48GB
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
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

plp = os.path.abspath('../../LiveOcean/plotting')
if plp not in sys.path:
    sys.path.append(plp)
import pfun

# set sample size
sampsize = 100

Ldir = Lfun.Lstart()

indir = '/data1/bbartos/LiveOcean_output/tracks/'
datadir = '/data1/bbartos/LiveOcean_data/tracker/'

# choose the run directory
print('\n%s\n' % '** Choose mooring file to plot **')
d_list_raw = os.listdir(indir)
d_list = []
for d in d_list_raw:
#    if d[-4:] == 'days':
        d_list.append(d)
Ndt = len(d_list)
for ndt in range(Ndt):
    print(str(ndt) + ': ' + d_list[ndt])
my_ndt = int(input('-- Input number -- '))
dirname = d_list[my_ndt] + '/'

# create the list of run files
m_list_raw = os.listdir(indir + dirname)
m_list = []
for m in m_list_raw:
    m_list.append(m)
Npt = len(m_list)
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
my_ndt = int(input('-- Input number (99 for all) -- '))
if my_ndt == 99:
    pass
else:
    m_list = [m_list[my_ndt],]

# output directory
odir00 = Ldir['parent'] + 'ptools_output/'
Lfun.make_dir(odir00)
odir0 = odir00 + 'tracks_bb/'
Lfun.make_dir(odir0)
outdir = odir0 + dirname
Lfun.make_dir(outdir)

# retrieve experimental data
exdf = pd.read_csv(datadir + 'rockfish_latlon.csv', index_col = 0)

# create plots for each run in the run directory
# if my_ndt == 99:
plt.ioff() # use this to supress plot output

for inname in m_list:
    
    # compile list of day files
    p_list = os.listdir(indir + dirname + inname)
    p_list.sort()
    # run through all days, concatenating the P dictionary in each
    counter = 0
    P = dict()
    for p in p_list:
        if counter == 0:
            # day 0 contains P, Ldir, and the grid data
            Pp, G, S, PLdir = pickle.load( open( indir + dirname + inname + '/' + p, 'rb' ) )
            # subsample particles from 100,000 particles
            samp = np.arange(5, Pp['lon'].shape[1], sampsize)
            for k in Pp.keys():
                if k in ['ot','age']:
                    P[k] = Pp[k]
                else:
                    P[k] = Pp[k][:,samp]
        else:
            # non-zero days only contain P and Ldir
            # first row overlaps with last row of previous day, so we remove it
            Pp, PLdir = pickle.load( open( indir + dirname + inname + '/' + p, 'rb' ) )
            for k in Pp.keys():
                if k in ['ot','age']:
                    P[k] = np.concatenate((P[k], Pp[k][1:]), axis=0)
                else:
                    P[k] = np.concatenate((P[k], Pp[k][1:,samp]), axis=0)
        counter += 1
        print('Finished ' + p)
        sys.stdout.flush()
    
    # set number of times and points
    NT, NP = P['lon'].shape
    # set map range
    aa = [np.min(P['lon'])-.1, np.max(P['lon'])+.1, np.min(P['lat'])-.1, 
          np.max(P['lat'])+.1]
    
    # get coastline
    cmat = matfun.loadmat(fn_coast)
    
    # retrieve relevant experimental information
    exind = inname[25:28]
    species = exdf['species'].loc[exind]
    location = exdf['Site'].loc[exind]
    
    # PLOTTING
    
    fig = plt.figure(figsize=(16,8))
    plt.suptitle(species+' rockfish released from '+location+' (2006)')
    ax = fig.add_subplot(1,2,1)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    
    # Coastline
    pfun.add_coast(ax)
    # Set axis limits
    ax.axis(aa)
    # Configure axis scales
    pfun.dar(ax)
    # Make lat/lon grid
    ax.grid()
    
    # tracks
    ax.plot(P['lon'][:], P['lat'][:], '-r', linewidth=0.5, alpha=0.5)
    # ending points
    ax.plot(P['lon'][-1,:],P['lat'][-1,:],'ob', markersize=4, label='End')
    # starting points
    ax.plot(P['lon'][0,:], P['lat'][0,:], '^y', markersize=10, label='Start')
    ax.legend()
        
    
    # TIME SERIES
    tdays = (P['ot'] - P['ot'][0])/86400.
    # Depth
    ax = fig.add_subplot(3,2,2)
    ax.plot(tdays, P['z'],'-', alpha=0.25)
    ax.set_ylabel('Z (m)')

    # Distance from Start
    dis = np.zeros(P['z'].shape)
    for tind in np.arange(1,len(tdays)):
    # change in lat/lon and total distance=sqrt(lat^2+lon^2)
        lat_dis = (P['lat'][tind,:] - P['lat'][tind-1,:]) * 111
        lon_dis = ((P['lon'][tind,:] - P['lon'][tind-1,:]) * 111 *
                        np.cos(np.nanmean(P['lat'][tind,:])*np.pi/180))
        dis[tind,:] = dis[tind-1,:] + np.sqrt(lat_dis**2 + lon_dis**2)
    # remove dead larvae
    dis_alive = np.where(np.isfinite(P['z']), dis, np.nan)
    ax = fig.add_subplot(3,2,4)
    ax.plot(tdays, dis_alive, '-', alpha=0.25)
    ax.set_ylabel('Distance Traveled (km)')
    ax.grid()
    
    # Net Distance from Start
    ndis = np.zeros(P['z'].shape)
    for tind in np.arange(1,len(tdays)):
    # change in lat/lon and total distance=sqrt(lat^2+lon^2)
        lat_ndis = (P['lat'][tind,:] - P['lat'][0,:]) * 111
        lon_ndis = ((P['lon'][tind,:] - P['lon'][0,:]) * 111 *
                        np.cos(np.nanmean(P['lat'][tind,:])*np.pi/180))
        ndis[tind,:] = np.sqrt(lat_ndis**2 + lon_ndis**2)
    # remove dead larvae
    ndis_alive = np.where(np.isfinite(P['z']), ndis, np.nan)
    ax = fig.add_subplot(3,2,6)
    ax.plot(tdays, ndis_alive, '-', alpha=0.25)
    ax.set_ylabel('Distance From Release(km)')
    ax.set_xlabel('Days')
    ax.grid()

    # save figures
    outfn = outdir + inname + '.png'
    plt.savefig(outfn)
    
    if my_ndt == 99:
        plt.close('all')
    else:
        plt.show()
        
plt.ion()

