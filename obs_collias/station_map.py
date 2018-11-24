#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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

dir0 = Ldir['parent'] + 'ptools_data/collias/'
# load processed station info and data
sta_df = pd.read_pickle(dir0 + 'sta_df.p')

# plotting
plt.close('all')
fig = plt.figure(figsize=(10,12))
ax = fig.add_subplot(111)

sta_df.plot(x='Longitude', y = 'Latitude', style='*r', ax=ax, legend=False)

for sta in sta_df.index:
    x = sta_df.loc[sta,'Longitude']
    y = sta_df.loc[sta,'Latitude']
    if isinstance(x,float) and isinstance(y,float):
        ax.text(x,y,sta, fontsize=6)
    else:
        print('')
        print(sta)
        print(x)
    
pfun.add_coast(ax)
pfun.dar(ax)

ax.set_title('Collias Stations')
ax.set_xlim(-125, -122)
ax.set_ylim(47, 49)

plt.show()