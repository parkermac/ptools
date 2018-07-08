#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 15:15:22 2018

@author: pm7

Plot locations of Stations.

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()

pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

import pandas as pd
import matplotlib.pyplot as plt

dir0 = Ldir['parent'] + 'ptools_data/ecology/'
# load processed station info and data
sta_df = pd.read_pickle(dir0 + 'sta_df.p')

# add Canadian data

dir1 = Ldir['parent'] + 'ptools_data/canada/'
# load processed station info and data
sta_df_ca = pd.read_pickle(dir1 + 'sta_df.p')

sta_df = pd.concat((sta_df, sta_df_ca))

# plotting
plt.close('all')
fig1 = plt.figure(figsize=(8,8))
fig2 = plt.figure(figsize=(8,8))
ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)

for station in sta_df.index:
    lon = sta_df.loc[station, 'Longitude']
    lat = sta_df.loc[station, 'Latitude']
    
    if sta_df.loc[station, 'Basin'] in ['Grays Harbor', 'Willapa Bay']:
        ax = ax1
    else:
        ax = ax2
    ax.plot(lon, lat, '*r')
    ax.text(lon+.01, lat, station, color='b', fontsize=8, fontweight='bold')
    
pfun.add_coast(ax1)
pfun.dar(ax1)

pfun.add_coast(ax2)
pfun.dar(ax2)

# Coastal Estuaries
ax1.set_title('Coastal Estuaries')
ax1.set_xlim(-124.3, -123.6)
ax1.set_ylim(46.3, 47.1)

# Puget Sound
ax2.set_title('Salish Sea')
ax2.set_xlim(-124, -122)
ax2.set_ylim(47, 49.5)

plt.show()

