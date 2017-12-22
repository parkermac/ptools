# -*- coding: utf-8 -*-
"""
Plot results of tracker using the new .nc files.
"""

# setup
import os
import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

plp = os.path.abspath('../../LiveOcean/plotting')
if plp not in sys.path:
    sys.path.append(plp)
import pfun

import numpy as np
import netCDF4 as nc
from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages


if Ldir['env'] == 'pm_mac': # mac version
    pass
elif Ldir['env'] == 'fjord': # fjord version
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt

# create the list of run files
indir = odir00 = Ldir['parent'] + 'ptools_output/rockfish/'
m_list_raw = os.listdir(indir)
m_list = []
for m in m_list_raw:
    if ('.nc' in m) and ('grid' not in m):
        m_list.append(m)
m_list.sort()
if False:
    Npt = len(m_list)
    for npt in range(Npt):
        print(str(npt) + ': ' + m_list[npt])
else:
    my_ndt = 0
inname = m_list[my_ndt]
ds = nc.Dataset(indir + inname)

# get data
NP = 100
lon = ds['lon'][:,:NP]
lat = ds['lat'][:,:NP]
z = ds['z'][:,:NP]
age = ds['age'][:,:NP]
ot = ds['ot'][:]
ds.close()

# PLOTTING
plt.close('all')
fig = plt.figure(figsize=(16,8))

# MAP
ax = fig.add_subplot(1,2,1)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
pfun.add_coast(ax)
aa = [-125, -122, 47, 49]
ax.axis(aa)
pfun.dar(ax)
ax.grid()

ax.plot(lon, lat, '-', linewidth=0.5, alpha=0.5)
# ending points
ax.plot(lon[-1], lat[-1],'ob', markersize=3, label='End')
# starting points
ax.plot(lon[0,0], lat[0,0], '*m', markersize=5, label='Start')
    
# TIME SERIES
tdays = (ot - ot[0])/86400.
ax = fig.add_subplot(1,2,2)
ax.plot(tdays, z,'-', alpha=0.25)
ax.set_ylabel('Z (m)')

# save or plot figures
outfn = indir + inname.strip('.nc') + '.png'
if Ldir['env'] == 'fjo':
    plt.savefig(outfn)
elif Ldir['env'] == 'pm_mac':
    plt.show()

# EXTRA

if False:
    
    # retrieve experimental data
    if Ldir['env'] == 'fjo':
        datadir = '/data1/bbartos/LiveOcean_data/tracker/'
    elif Ldir['env'] == 'pm_mac':
        datadir = Ldir['parent'] + 'ptools_data/rockfish/'
    exdf = pd.read_csv(datadir + 'rockfish_latlon.csv', index_col = 0)
    exind = inname[-6:-3]
    species = exdf['species'].loc[exind]
    location = exdf['Site'].loc[exind]
    

