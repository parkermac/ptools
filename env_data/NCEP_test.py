# -*- coding: utf-8 -*-
"""
Code to test NCEP sunshine data.

"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta

pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart('cascadia1','base')

dt0 = datetime(1965,1,1)
dt1 = datetime(1984,1,1)

#%% NCEP Sunshine Re-analaysis
    
dir0 = '/Users/PM5/Documents/'
indir0 = dir0 + 'ptools_output/env_data/'
ncep_indir = indir0 + 'reanalysis/'
ncep_fn = 'df_NCEP_dswrf_Seattle.p'

ncep_df = pd.read_pickle(ncep_indir + ncep_fn)

ncep_df = ncep_df[(ncep_df.index>=dt0) & (ncep_df.index<=dt1)]

#%% Weather data

dir0 = '/Users/PM5/Documents/tools_data/obs_data/met/'
wea_fn = dir0 + 'seatac_1948-2015_v2.csv'
# read csv data into a data frame
  
df = pd.read_csv(wea_fn, index_col = 'DATE', parse_dates = ['DATE'])

df[df == -9999] = np.nan

cols = list(df.columns)
sta_name = df.ix[0,'STATION_NAME']
sta_id = df.ix[0,'STATION']

# make a list of tuples with:
# (variable name, scaling factor, label text, low, high)

vtup_list = [
    ('TMAX',0.1,'Tmax (degC)',0,30),
    ('TMIN',0.1,'Tmin (degC)',-10,20),
    ('PRCP',0.1,'Precip (mm/day)',0,10),
    ('ACMH',1.0,'Clouds (obs) %',0,100),
    #('ACSH',1.0,'Daytime Clouds (obs) %',0,100),
    ('PSUN',1.0,'% of Possible Sunshine',0,100),
    ('TSUN',1/60.,'Daily Sunshine (hours)',0,12),
    #('AWND',0.1,'Average Daily Windspeed (m/s)',0,5),
    #('WDF2',1.0,'Dir of 2-minute Wind (deg)',0,360)
    ]

vname_list = []
for tup in vtup_list:
    vname_list.append(tup[0])

df = df[df.index<datetime(2015,1,1)]

df = df[vname_list]

for tup in vtup_list:
    vn = tup[0]
    df[vn] = df[vn]*tup[1]

# get monthly means
# the period 'MS' does monthly with the timestamp at the START of the month
# so in theory this would have a 15 day lead...
dfm = df.resample('MS', how='mean',
                  closed = 'left', label='left')

dfm = dfm[(dfm.index>=dt0) & (dfm.index<=dt1)]

# climatology
# mo_fact is to try to account for the fact that there is more
# daylight in the summer
mo_fact = 1. + 0.5*np.cos(np.linspace(-np.pi,np.pi,12))
dfma = dfm.copy()
for mm in range(1,13):
    im = dfm.index.month == mm
    dfma[im] = mo_fact[mm-1] * (dfm[im] - dfm[im].mean())
    
#%% Combining

anoms = pd.DataFrame(index=ncep_df.index)
anoms['ncep'] = ncep_df['anom']
anoms['psun'] = dfma['PSUN']
anoms['tsun'] = dfma['TSUN']
    
#%% Plotting

plt.close('all')

fig_size = (25,10)
fig = plt.figure(figsize=fig_size)

ax = plt.subplot2grid((1,3), (0,0), colspan=2)                        
ax.plot(anoms['ncep'], '-b')
ax.plot(anoms['psun'], '-r')
ax.grid()
ax.set_xlabel('Date')
ax.text(.1, .9, 'NCEP Downward SW Rad. Anomaly ($W$ $m^{2}$)',
        fontsize=16, color='b', transform=ax.transAxes)
ax.text(.1, .1, 'SeaTac Sunshine Anomaly (%)',
        fontsize=16, color='r', transform=ax.transAxes)


ax = plt.subplot2grid((1,3), (0,2), colspan=1) 
ax.plot(anoms['ncep'],anoms['psun'], '*k')
ax.grid()
ax.set_xlabel('SeaTac Sunshine Anomaly (%)')
ax.set_ylabel('NCEP Downward SW Rad. Anomaly ($W$ $m^{2}$)')


plt.show()

