"""
Code to process Hakai casts.

The two files each correspond to a single station, the first a N SoG,
and the second at the NW end of Johnstone Strait.
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
from datetime import datetime, timedelta
import seawater as sw

indir = '../../ptools_data/hakai/'
fn1 = indir + '/raw/ctd-bulk-1539111774185.csv'
fn2 = indir + '/raw/ctd-bulk-1539112518591.csv'

a1 = pd.read_csv(fn1)
a2 = pd.read_csv(fn2)

# print info about station location
print('a1: lon = %0.1f (+/- %0.5f) lat = %0.1f (+/- %0.5f)' % 
    (a1['Longitude'].mean(),a1['Longitude'].std(),
    a1['Latitude'].mean(),a1['Latitude'].std()))
print('a2: lon = %0.1f (+/- %0.5f) lat = %0.1f (+/- %0.5f)' % 
    (a2['Longitude'].mean(),a2['Longitude'].std(),
    a2['Latitude'].mean(),a2['Latitude'].std()))
    
# get cast time info
cn_list = set(a1['Cast PK'])
c1 = pd.DataFrame()
for cn in cn_list:
    c = a1[a1['Cast PK']==cn]
    c1.loc[cn,'Date'] = c['Start time'].values[0]
    
cn_list = set(a2['Cast PK'])
c2 = pd.DataFrame()
for cn in cn_list:
    c = a2[a2['Cast PK']==cn]
    c2.loc[cn,'Date'] = c['Start time'].values[0]
    
c1 = c1.sort_values('Date')
c2 = c2.sort_values('Date')

# PLOTTING
plt.close('all')

fig = plt.figure(figsize=(16,10))

# map
ax = fig.add_subplot(211)
a1.plot(x='Longitude', y='Latitude', style='.r', ax = ax, legend=False)
a2.plot(x='Longitude', y='Latitude', style='.b', ax = ax, legend=False)
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis([-127.5, -124.5, 50, 50.8])
ax.set_title('Hakai Stations')
ax.set_xlabel('Longitude (deg)')
ax.set_ylabel('Latitude (deg)')

# cast data
axs = fig.add_subplot(245) # salinity
axt = fig.add_subplot(246) # temperature
axo = fig.add_subplot(247) # oxygen
axd = fig.add_subplot(248) # density

cn_list = set(a1['Cast PK'])
for cn in cn_list:
    c = a1[a1['Cast PK']==cn].copy()
    c.plot(x='Salinity (PSU)', y='Depth (m)',style='r-', ax=axs, legend=False)
    c.plot(x='Temperature (deg C)', y='Depth (m)',style='r-', ax=axt, legend=False)
    o = (1.42903 * 1000 / 32) * c['Dissolved O2 (mL/L)'].values
    c.loc[:,'DO (uM)'] = o
    c.plot(x='DO (uM)', y='Depth (m)',style='r-', ax=axo, legend=False)
    p = c['Pressure (dbar)'].values
    s = c['Salinity (PSU)'].values
    t = c['Temperature (deg C)'].values
    pden = sw.pden(s, t, p, pr=0)
    sig = pden-1000
    c.loc[:,'Sigma (kg/m3)'] = sig
    c.plot(x='Sigma (kg/m3)', y='Depth (m)',style='r-', ax=axd, legend=False)

cn_list = set(a2['Cast PK'])
for cn in cn_list:
    c = a2[a2['Cast PK']==cn].copy()
    c.plot(x='Salinity (PSU)', y='Depth (m)',style='b-', ax=axs, legend=False)
    c.plot(x='Temperature (deg C)', y='Depth (m)',style='b-', ax=axt, legend=False)
    o = (1.42903 * 1000 / 32) * c['Dissolved O2 (mL/L)'].values
    c.loc[:,'DO (uM)'] = o
    c.plot(x='DO (uM)', y='Depth (m)',style='b-', ax=axo, legend=False)
    p = c['Pressure (dbar)'].values
    s = c['Salinity (PSU)'].values
    t = c['Temperature (deg C)'].values
    pden = sw.pden(s, t, p, pr=0)
    sig = pden-1000
    c.loc[:,'Sigma (kg/m3)'] = sig
    c.plot(x='Sigma (kg/m3)', y='Depth (m)',style='b-', ax=axd, legend=False)

axs.grid(True)
axt.grid(True)
axo.grid(True)
axd.grid(True)

dmax = 200
axs.set_ylim(dmax,0)
axt.set_ylim(dmax,0)
axo.set_ylim(dmax,0)
axd.set_ylim(dmax,0)

plt.show()
    
    
    