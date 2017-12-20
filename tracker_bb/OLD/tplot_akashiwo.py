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
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
my_npt = int(input('-- Input number -- '))
inname = m_list[my_npt]

# output directory
outdir = indir + 'plots/' + dirname
Lfun.make_dir(outdir)
#Lfun.make_dir(outdir, clean=True) # use this to clear previous plots

# supress plot output
plt.ioff()

# load run dictionaries
Pfull, G, S, PLdir = pickle.load( open( indir + dirname + inname, 'rb' ) )

NT, NP = Pfull['lon'].shape

# set map limits
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

# number of times (snapshots) to plot
ntotal = 4

bm_dict = dict()
P_dict = dict()

# set index to plot
index_list = []
for frac in np.arange(1,ntotal+1):
    index_num = int(NT*frac/ntotal)
    index_list.append(index_num)

# plot each snapshot separately
for ind in index_list:    
    
    # set date
    oceant = Pfull['ot'][ind-1]
    full_date = str(datetime(2009, 1, 1) + timedelta(seconds=oceant))
    date = full_date[:10]
    
    # define subdirectory of Pfull
    P = dict()
    for key in Pfull:
        P[key] = Pfull[key][0:ind]
    
    P_dict[date] = P
    
    # create figure
    fig = plt.figure(figsize=(16,16))
    ax = plt.gca()
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(inname + '_' + date)
    
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
    
    # plot all tracks
    for ii in range(NP):
        ax.plot(P['lon'][0, ii],P['lat'][0, ii],'ok', alpha = .1)
    # mask of tracks ending at the shore
    cc = 0
    beach_mask = []
    for j in range(P['u'].shape[1]):
        if P['u'][-1,j] == 0 and P['v'][-1,j] == 0:
            beach_mask.append(True)
        else:
            beach_mask.append(False)
    beach_mask = np.array(beach_mask)
#    beach_mask = P['u'][-1,:] == 0 and P['v'][-1,:] == 0

    bm_dict[date] = beach_mask
    
    print(str(sum(beach_mask)) + ' tracks on shore by ' + date)
    # start points
    ax.plot(P['lon'][0,beach_mask], P['lat'][0,beach_mask], 'ob', markersize=5, label='Start')
    # end points
    ax.plot(P['lon'][-1,beach_mask], P['lat'][-1,beach_mask], 'or', markersize=5, label='End')
    print(str(len(P['lon'][-1,beach_mask])) + ' endpoints')
    ax.legend()
    
    # contour legend
    ax.text(.85, .25, 'Depth Contours', horizontalalignment='center', transform=ax.transAxes, color='g')
    dd = 1
    for d in depth_levs:
        ax.text(.85, .25 - .03*dd, str(d), horizontalalignment='center', transform=ax.transAxes, color='g')
        dd += 1
    
    # save figures
    outfn = outdir + inname[:-2] + '_' + date + '.png'
    plt.savefig(outfn)
    
# turn output back on
plt.close('all')
plt.ion()