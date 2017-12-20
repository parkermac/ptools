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
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

plp = os.path.abspath('../plotting')
if plp not in sys.path:
    sys.path.append(plp)
#import pfun

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'tracks/'
fn_coast = Ldir['data'] + 'coast/pnw_coast_combined.mat'

# choose the type of plot to make
print('\n%s\n' % '** Choose mooring file to plot **')
m_list_raw = os.listdir(indir)
m_list = []
for m in m_list_raw:
    if m[-2:] == '.p':
        m_list.append(m)
Npt = len(m_list)
m_dict = dict(zip(range(Npt), m_list))
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
my_npt = int(input('-- Input number -- '))
inname = m_dict[my_npt]

import pickle  # python 3
P, G, S, PLdir = pickle.load( open( indir + inname, 'rb' ) )

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
else:
    aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]
depth_levs = [100, 200, 500, 1000, 2000, 3000]

# PLOTTING
#plt.close('all')
if 'akashiwo' in inname:
    fig = plt.figure()
    ax = plt.gca()
else:
    fig = plt.figure(figsize=(16,8))
    ax = fig.add_subplot(1,2,1)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title(inname)

# construct basemap
m = Basemap(llcrnrlon=int(aa[0])-.5, llcrnrlat=int(aa[2])-.25, urcrnrlon=int(aa[1]), 
            urcrnrlat=int(aa[3]), projection='cyl', resolution='i')
m.drawcoastlines()
m.fillcontinents(color='grey', alpha=0.25)
m.drawcountries()
m.drawstates()
m.drawmapboundary()
m.drawmeridians(np.arange(int(aa[0]),int(aa[1]),2), labels=[0,0,0,1])
m.drawparallels(np.arange(int(aa[2]),int(aa[3]),2), labels=[1,0,0,0])

# density contours
m.contour(G['lon_rho'], G['lat_rho'], G['h'], depth_levs, colors='k')

# plot tracks
cs_divider = -.4
csd_text = str(int(np.abs(100.*cs_divider)))
ii = 0
for cs in P['cs'][NT-1,:]:
    if cs < cs_divider:
        m.plot(P['lon'][:, ii], P['lat'][:, ii], '-b', alpha = .4, latlon=True)
    else:
        m.plot(P['lon'][:, ii], P['lat'][:, ii], '-r', alpha = .4, latlon=True)
    ii += 1

# starting points
ii = 0
for cs in P['cs'][NT-1,:]:
    if cs < cs_divider:
        m.plot(P['lon'][0,ii], P['lat'][0,ii], 'ob', markersize=5, alpha = .4, 
               markeredgecolor='b', latlon=True)
    else:
        m.plot(P['lon'][0,ii], P['lat'][0,ii], 'or', markersize=5, alpha = .4, 
               markeredgecolor='r', latlon=True)
    ii += 1

# ending points
#m.plot(P['lon'][-1,:], P['lat'][-1,:], 'yo', markersize=5, latlon=True)

# contour legend
ax.text(.85, .25, 'Depth Contours', horizontalalignment='center', 
         transform=ax.transAxes, color='k')
dd = 1
for d in depth_levs:
    ax.text(.85, .25 - .03*dd, str(d), horizontalalignment='center', 
           transform=ax.transAxes, color='k')
    dd += 1

# text for depth percentages:
if PLdir['surface'] == True:
    pass
else:
    ax.text(.95,.9, 'End depth above '+csd_text+'%', horizontalalignment='right', 
        transform=ax.transAxes, color='r', fontsize=16)
    ax.text(.95,.8, 'End depth below '+csd_text+'%', horizontalalignment='right', 
        transform=ax.transAxes, color='b', fontsize=16)

# to remove time series use:
if 'akashiwo' in inname:
    pass
else:
    # TIME SERIES
    tdays = (P['ot'] - P['ot'][0])/86400.
    
    # u velocity
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

    # v velocity
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

    # depth
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
    plt.show()