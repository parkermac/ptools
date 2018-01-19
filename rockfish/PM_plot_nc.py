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
import pickle

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

if Ldir['lo_env'] == 'pm_mac': # mac version
    pass
elif Ldir['lo_env'] == 'pm_fjord': # fjord version
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt

# create the list of run files
indir = odir00 = Ldir['parent'] + 'ptools_output/rockfish/'
ex_list_raw = os.listdir(indir)
ex_list = []
for ex in ex_list_raw:
    if ('.nc' in ex) and ('grid' not in ex):
        ex_list.append(ex)
ex_list.sort()

# testing
ex = 'rockfish_2006_Experiment_4_4.nc'

# get data
ds = nc.Dataset(indir + ex)
# tracks are stored (time, particle)
Lon = ds['lon'][:]
Lat = ds['lat'][:]
Z = ds['z'][:]
H = ds['h'][:]
Age = ds['age'][:]
ot = ds['ot'][:]
ds.close()

# subsample the number of particles
pstep = 100
lon = Lon[:,::pstep]
lat = Lat[:,::pstep]
z = Z[:,::pstep]
h = H[:,::pstep]
age = Age[:,::pstep]
NT, NP = lon.shape

# find time indices that define an
# age range (e.g. 0-120 days)
ndays0 = 0
ndays1 = 120
age_ind0 = (age == ndays0).sum(axis=0)
age_ind1 = (age <= ndays1).sum(axis=0)

# repackage to all have the same age
nt = age_ind1[0] - age_ind0[0]
llon = np.nan * np.ones((nt, NP))
llat = np.nan * np.ones((nt, NP))
for np in range(NP):
    nt0 = age_ind1[0]
    nt1 = nt0 + nt
    llon[:nt,np] = lon[nt0:nt1, np]
    llat[:nt,np] = lat[nt0:nt1, np]

# PLOTTING
plt.close('all')
fig = plt.figure(figsize=(13,7))

# MAP
ax = fig.add_subplot(1,2,1)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
pfun.add_coast(ax)
# automatically plot region of particles, with padding
pad = .02
aa = [lon.min() - pad, lon.max() + pad,
    lat.min() - pad, lat.max() + pad]
#aa = [-123.5, -122, 47, 48.5]
ax.axis(aa)
pfun.dar(ax)
ax.grid()
ax.plot(llon, llat,
    '-', linewidth=0.5, alpha=0.5)
# ax.plot(lon[age_ind0[np],np], lat[age_ind0[np],np],
#     '*m', markersize=5) # start
# ax.plot(lon[age_ind1[np],np], lat[age_ind1[np],np],
#     'ob', markersize=3) # end

# TIME SERIES
tdays = (ot - ot[0])/86400.
ax = fig.add_subplot(2,2,2)
ax.plot(tdays, z,'-', alpha=0.25)
ax.set_ylabel('Z (m)')

# ax = fig.add_subplot(2,2,4)
# ax.plot(age, z,'.')

# save or plot figures
if Ldir['lo_env'] == 'pm_fjord':
    out_fn = indir + ex.strip('.nc') + '.png'
    plt.savefig(out_fn)
elif Ldir['lo_env'] == 'pm_mac':
    plt.show()

if False:
    # retrieve experimental data
    if Ldir['lo_env'] == 'pm_fjord':
        datadir = '/data1/bbartos/LiveOcean_data/tracker/'
    elif Ldir['lo_env'] == 'pm_mac':
        datadir = Ldir['parent'] + 'ptools_data/rockfish/'
    exdf = pd.read_csv(datadir + 'rockfish_latlon.csv', index_col = 0)
    exind = inname[-6:-3]
    species = exdf['species'].loc[exind]
    location = exdf['Site'].loc[exind]
    

