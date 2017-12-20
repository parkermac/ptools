# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:03:12 2016

@author: PM5
"""

"""
Plot results of tracker.
"""

# setup
import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import matfun
import pickle
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

plp = os.path.abspath('../plotting')
if plp not in sys.path:
    sys.path.append(plp)
import pfun

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'tracks/'
fn_coast = Ldir['data'] + 'coast/pnw_coast_combined.mat'

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
outdir = indir + 'plots/' + dirname
Lfun.make_dir(outdir)
#Lfun.make_dir(outdir, clean=True) # use this to clear previous plots

# create plots for each run in the run directory
if my_ndt == 99:
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
            for k in Pp.keys():
                P[k] = Pp[k]
        else:
            # non-zero days only contain P and Ldir
            # first row overlaps with last row of previous day, so we remove it
            Pp, PLdir = pickle.load( open( indir + dirname + inname + '/' + p, 'rb' ) )
            for k in Pp.keys():
                if k == 'ot':
                    P[k] = np.concatenate((P[k], Pp[k][1:]), axis=0)
                else:
                    P[k] = np.concatenate((P[k], Pp[k][1:,:]), axis=0)
        counter += 1    
    
    NT, NP = P['lon'].shape
    
    lonp = G['lon_psi']
    latp = G['lat_psi']
    if 'jdf' in inname:
        aa = [-126, -123.5, 47, 49]
    elif 'cr' in inname:
        aa = [-125, -122.5, 45, 48]
    elif 'test' in inname:
        #aa = [-126.5, -125, 45, 46.2]
        aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]
    elif 'akashiwo' in inname:
        aa = [lonp.min(), lonp.max(), latp.min()-0.25, latp.max()]
    else:
        aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]
    depth_levs = [100, 200, 500, 1000, 2000, 3000]
    
    # get coastline
    cmat = matfun.loadmat(fn_coast)
    
    # PLOTTING
    
    #plt.close()
    
    # to not plot time series, add ic_name to this and the same structure below:
    if 'akashiwo' in inname:
        fig = plt.figure(figsize=(16,16))
        ax = plt.gca()
    else:
        fig = plt.figure(figsize=(16,8))
        ax = fig.add_subplot(1,2,1)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(inname)
    
    
    # MAP OF TRACKS
    
    # Depth Contours
    ax.contour(G['lon_rho'], G['lat_rho'], G['h'], depth_levs, colors='g')
    # Coastline
    ax.plot(cmat['lon'], cmat['lat'], '-k', linewidth=.5)
    # Set axis limits
    ax.axis(aa)
    # Configure axis scales
    pfun.dar(ax)
    # Make lat/lon grid
    ax.grid()
    
    # depth separation
    cs_divider = -.4
    csd_text = str(int(np.abs(100.*cs_divider)))
    
    # plotting tracks
    # to highlight tracks ending at the shore, add ic_name:
    if 'akashiwo' in inname:
        # plot all tracks
        for ii in range(NP):
            ax.plot(P['lon'][0, ii],P['lat'][0, ii],'ok', alpha = .1)
        # mask of tracks ending at the shore
        beach_mask = P['u'][-1,:] == 0
        # tracks
        ax.plot(P['lon'][:,beach_mask], P['lat'][:,beach_mask], '-r', linewidth=1)
        # start points
        ax.plot(P['lon'][0,beach_mask], P['lat'][0,beach_mask], 'ob', markersize=5, label='Start')
        # end points
        ax.plot(P['lon'][-1,beach_mask], P['lat'][-1,beach_mask], 'or', markersize=5, label='End')
        ax.legend()
    else:
        ii = 0
        for cs in P['cs'][NT-1,:]:
            if cs < cs_divider:
                ax.plot(P['lon'][:, ii],P['lat'][:, ii],'-b', alpha = .4)
            else:
                ax.plot(P['lon'][:, ii],P['lat'][:, ii],'-r', alpha = .4)
            ii += 1
        # starting points
        ii = 0
        for cs in P['cs'][NT-1,:]:
            if cs < cs_divider:
                ax.plot(P['lon'][0,ii], P['lat'][0,ii], 'ob', markersize=5, alpha = .4, markeredgecolor='b')
            else:
                ax.plot(P['lon'][0,ii], P['lat'][0,ii], 'or', markersize=5, alpha = .4, markeredgecolor='r')
            ii += 1
    
        # ending points
        #ax.plot(P['lon'][-1,:],P['lat'][-1,:],'y*',markersize=20)
    
    # contour legend
    ax.text(.85, .25, 'Depth Contours', horizontalalignment='center', transform=ax.transAxes, color='g')
    dd = 1
    for d in depth_levs:
        ax.text(.85, .25 - .03*dd, str(d), horizontalalignment='center', transform=ax.transAxes, color='g')
        dd += 1
    
    # text for depth percentages:
    if PLdir['surface'] == True:
        pass
    else:
        ax.text(.95, .9, 'End depth above '+csd_text+'%', horizontalalignment='right', 
            transform=ax.transAxes, color='r', fontsize=16)
        ax.text(.95, .8, 'End depth below '+csd_text+'%', horizontalalignment='right', 
            transform=ax.transAxes, color='b', fontsize=16)
    
    # to not plot time series, add ic_name to this and the same structure above:
    if 'akashiwo' in inname:
        pass
    
    else:
        # TIME SERIES
        tdays = (P['ot'] - P['ot'][0])/86400.
    
        ax = fig.add_subplot(3,2,2)
        ii = 0
        for cs in P['cs'][NT-1,:]:
            if cs < cs_divider:
                ax.plot(tdays, P['u'][:,ii],'-b')
            else:
                ax.plot(tdays, P['u'][:,ii],'-r')
            ii += 1
        ax.set_ylabel('U $m s^{-1}$')
        ax.set_ylim(-.8, .8)
        ax.grid()
    
        ax = fig.add_subplot(3,2,4)
        ii = 0
        for cs in P['cs'][NT-1,:]:
            if cs < cs_divider:
                ax.plot(tdays, P['v'][:,ii],'-b')
            else:
                ax.plot(tdays, P['v'][:,ii],'-r')
            ii += 1
        ax.set_ylabel('V $m s^{-1}$')
        ax.set_ylim(-.8, .8)
        ax.grid()
    
        ax = fig.add_subplot(3,2,6)
        ii = 0
        for cs in P['cs'][NT-1,:]:
            if cs < cs_divider:
                ax.plot(tdays, P['z'][:,ii],'-b')
            else:
                ax.plot(tdays, P['z'][:,ii],'-r')
            ii += 1
        ax.set_xlabel('Days')
        ax.set_ylabel('Z (m)')
    
    # save figures
    outfn = outdir + inname[:-2] + '.png'
    plt.savefig(outfn)
    
    if my_ndt != 99:
        plt.show()

if my_ndt == 99:
    plt.close('all')
    plt.ion()
    