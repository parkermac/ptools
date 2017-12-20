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
Lfun.make_dir(indir)
fn_coast = Ldir['data'] + 'coast/pnw_coast_combined.mat'

# choose the run directory
print('\n%s\n' % '** Choose mooring file to plot **')
d_list_raw = os.listdir(indir)
d_list = []
for d in d_list_raw:
    if d[-4:] == 'days':
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
    if m[-2:] == '.p':
        m_list.append(m)
Npt = len(m_list)

# output directory
outdir0 = indir + 'plots/'
Lfun.make_dir(outdir0)
outdir = outdir0 + dirname
Lfun.make_dir(outdir)
#Lfun.make_dir(outdir, clean=True) # use this to clear previous plots

# create plots for each run in the run directory
#plt.ioff() # use this to supress plot output
for inname in m_list:
    
    P, G, S, PLdir = pickle.load( open( indir + dirname + inname, 'rb' ) )
    
    NT, NP = P['lon'].shape
    
    lonp = G['lon_psi']
    latp = G['lat_psi']
    aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]
    depth_levs = [100, 200, 500, 1000, 2000, 3000]
    
    # get coastline
    cmat = matfun.loadmat(fn_coast)
    
    # PLOTTING
    
    #plt.close()
    
    # to not plot time series, add ic_name to this and the same structure below:
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

    # TIME SERIES
    tdays = (P['ot'] - P['ot'][0])/86400.

    ax = fig.add_subplot(1,2,2)
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
#    plt.close()
    
#plt.ion()