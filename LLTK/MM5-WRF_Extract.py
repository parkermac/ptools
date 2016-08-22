# -*- coding: utf-8 -*-
"""
Code to extract and explore wind time series, for the LLTK project.

This is using the MM5-WRF archive I created as part of the RISE and PNWTOX
projects.

"""

# setup
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import h5py

pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart('cascadia1','base')
fn_coast = Ldir['data'] + 'coast/pnw_coast_combined.mat'
import zfun
import matfun

#%% naming

dir0 = '/Users/PM5/Documents/'
indir = dir0 + 'tools_data/forcing_data/atm/mm5/tinterp_three_hour/'
outdir = dir0 + 'ptools_output/env_data/mm5_wrf/'
try:
    os.mkdir(outdir)
except OSError:
    pass # assume OSError was raised because directory already exists

vn_list = ['u10r', 'v10r', 'swdown', 't2']
c_list = ['b', 'r'] # colors for lines
iLon_target = -123.
iLat_target = 48.3

fld_dict = dict()
ifld_dict = dict()

# first get a single spatial map
yr = 2005
mo = 6
fld_dict = dict()
# packed x,y,t!
for vn in vn_list:
    fn = indir + str(yr) + '/' + vn + '_' + str(mo) + '.mat'
    f = h5py.File(fn)   
    ff = dict()
    for item in f.keys():
        ff[item] = f[item][:]
        # print(item + ' ' + str(ff[item].shape))
    f.close()
    lon = ff['LON']
    lat = ff['LAT']                
    fld = ff['VAR'][:,:,0]    
    fld_dict[vn] = fld
Lon = lon[:,0]
Lat = lat[0,:]
ix = np.digitize(iLon_target, Lon)
iy = np.digitize(iLat_target, Lat)
iLon = Lon[ix]
iLat = Lat[iy]

# Then get a time series at the chosen location
ifld_dict = dict()

for vn in vn_list:
    ifld = np.array([])
    it_datetime = []
    for yr in range(2002,2010):    
        for mo in range(1,13):      
            fn = indir + str(yr) + '/' + vn + '_' + str(mo) + '.mat'
            f = h5py.File(fn)   
            this_ifld = f['VAR'][ix, iy, :]               
            this_it_datenum = f['TD'][:].squeeze()
            f.close()
            ifld = np.concatenate((ifld, this_ifld))
            # convert matlab datenum to python datetime
            for dn in this_it_datenum:
                it_datetime.append( datetime.fromordinal(dn.astype(int)) +
                    timedelta(days = np.mod(dn,1)) - timedelta(days = 366) )         
    ifld_dict[vn] = ifld

#%% pack results in a DataFrame

wind_df = pd.DataFrame(index=it_datetime)
for vn in vn_list:
    wind_df[vn] = ifld_dict[vn]
    
wind_df.to_pickle(outdir + 'mm5_wrf_jdfE.p')        

#%% Plotting

plt.close()

cmat = matfun.loadmat(fn_coast)

fig_size = (20, 10)
fig = plt.figure(figsize=fig_size)

counter = 0

for vn in vn_list:
    
    ax = plt.subplot2grid((3,len(vn_list)), (0,counter), rowspan=2)    
    axh = ax.pcolormesh(lon,lat,fld_dict[vn])
    plt.colorbar(axh)
    ax.plot(iLon, iLat, '*m', markersize=15)
    ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
    ax.set_xlim(Lon[0], Lon[-1])
    ax.set_ylim(Lat[0], Lat[-1])
    ax.set_title(vn)
    zfun.dar(ax)
    
    counter += 1
    
ax = plt.subplot2grid((3,1), (2,0), colspan=1)
wind_df['u10r'].plot(ax=ax)
wind_df['v10r'].plot(ax=ax)
ax.set_xlim(it_datetime[0], it_datetime[-1])   
plt.show()

