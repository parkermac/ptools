"""
Calculates validation statistics and other quantities for temperature,
salinity and DO from Ecology and EC CTD casts.  Uses all depths in 1 m
increments.

Good performance - under 30 sec per year.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart()
import zfun

# ***** User Edits

# +++ load Ecology and Canadian CTD cast data +++
dir0 = Ldir['parent'] + 'ptools_data/ecology/'
dir1 = Ldir['parent'] + 'ptools_data/canada/'
sta_df = pd.read_pickle(dir0 + 'sta_df.p')
sta_df_ca = pd.read_pickle(dir1 + 'sta_df.p')
sta_df = pd.concat((sta_df, sta_df_ca), sort=False)

Ldir['gtagex'] = 'cas6_v3_lo8b'

testing = False

# ***** End User Edits

TSO_all = pd.DataFrame()
TSO_dict = dict()

for year in [2017, 2018, 2019]:
    
    # input location
    in_dir = Ldir['parent'] + 'ptools_output/ecology/'
    in_fn = 'TSO_1m_casts_' + Ldir['gtagex'] + '_' + str(year) + '.p'
    
    TSO = pd.read_pickle(in_dir + in_fn)
    TSO_dict[year] = TSO
    TSO_all = pd.concat((TSO_all, TSO), sort=False)
    
    # statistics by year
    DS = TSO['Mod Salinity'] - TSO['Salinity']
    DT = TSO['Mod Temp. (deg C)'] - TSO['Temp. (deg C)']
    DO = TSO['Mod DO (mg L-1)'] - TSO['DO (mg L-1)']

    print('\n*** ' + str(year) + ' ***')
    print('Salinity Bias = %0.2f, RMSE = %0.2f' % (DS.mean(), np.sqrt((DS**2).mean())) )
    print('Temp. (deg C) Bias = %0.2f, RMSE = %0.2f' % (DT.mean(), np.sqrt((DT**2).mean())) )
    print('DO (mg L-1) Bias = %0.2f, RMSE = %0.2f' % (DO.mean(), np.sqrt((DO**2).mean())) )

# statistics all years
DS = TSO_all['Mod Salinity'] - TSO_all['Salinity']
DT = TSO_all['Mod Temp. (deg C)'] - TSO_all['Temp. (deg C)']
DO = TSO_all['Mod DO (mg L-1)'] - TSO_all['DO (mg L-1)']

print('\n*** ALL YEARS ***')
print('Salinity Bias = %0.2f, RMSE = %0.2f' % (DS.mean(), np.sqrt((DS**2).mean())) )
print('Temp. (deg C) Bias = %0.2f, RMSE = %0.2f' % (DT.mean(), np.sqrt((DT**2).mean())) )
print('DO (mg L-1) Bias = %0.2f, RMSE = %0.2f' % (DO.mean(), np.sqrt((DO**2).mean())) )

ts_df = pd.DataFrame(index=['SJF001','ADM002','ADM001','ADM003','PSB003','EAP001','CMB003',
        'GOR001','NSQ002','DNA001'],
        columns=['Longitude','Latitude', 'Distance (km)','S top','Mod S top' ,'S bot','Mod S bot',
        'T top','Mod T top', 'T bot','Mod T bot'])

x_list = []; y_list = []
for sn in ts_df.index:
        lon = sta_df.loc[sn,'Longitude']
        lat = sta_df.loc[sn,'Latitude']
        x, y = zfun.ll2xy(lon, lat, -125, 48)
        x_list.append(x)
        y_list.append(y)
xvec = np.array(x_list)/1000
yvec = np.array(y_list)/1000
xx = np.zeros_like(xvec)
yy = np.zeros_like(yvec)
xx[1:] = np.diff(xvec)
yy[1:] = np.diff(yvec)
dd = np.sqrt(xx**2 + yy**2)
dd = np.cumsum(dd)

# PLOTTING
plt.close('all')
fig = plt.figure(figsize=(14,11))

counter = 0
for year in [2017, 2018, 2019]:
    
    this_TSO = TSO_dict[year]

    ts_df['Distance (km)'] = dd
    for sn in ts_df.index:
        ts_df.loc[sn,'Longitude'] = sta_df.loc[sn,'Longitude']
        ts_df.loc[sn,'Latitude'] = sta_df.loc[sn,'Latitude']
        ts = this_TSO[this_TSO['Station']==sn]
        ts = ts[['Z','Salinity', 'Mod Salinity','Temp. (deg C)','Mod Temp. (deg C)']]
        ts = ts.dropna()
        ts_top = ts[ts['Z']>= -15.].mean()
        ts_bot = ts[ts['Z']< -15.].mean()
    
        ts_df.loc[sn,'S top'] = ts_top['Salinity']
        ts_df.loc[sn,'Mod S top'] = ts_top['Mod Salinity']
        ts_df.loc[sn,'S bot'] = ts_bot['Salinity']
        ts_df.loc[sn,'Mod S bot'] = ts_bot['Mod Salinity']
    
        ts_df.loc[sn,'T top'] = ts_top['Temp. (deg C)']
        ts_df.loc[sn,'Mod T top'] = ts_top['Mod Temp. (deg C)']
        ts_df.loc[sn,'T bot'] = ts_bot['Temp. (deg C)']
        ts_df.loc[sn,'Mod T bot'] = ts_bot['Mod Temp. (deg C)']

    ax = fig.add_subplot(3,2,counter*2 + 1)
    ts_df.plot(x='Distance (km)', y=['S top','Mod S top'], ax=ax, style=['-','--'], color='r')
    ts_df.plot(x='Distance (km)', y=['S bot','Mod S bot'], ax=ax, style=['-','--'], color='b')
    ax.set_ylim(28,32)
    
    ax = fig.add_subplot(3,2,counter*2 + 2)
    ts_df.plot(x='Distance (km)', y=['T top','Mod T top'], ax=ax, style=['-','--'], color='r')
    ts_df.plot(x='Distance (km)', y=['T bot','Mod T bot'], ax=ax, style=['-','--'], color='b')
    ax.set_ylim(8,12)
    
    counter += 1

plt.show()
