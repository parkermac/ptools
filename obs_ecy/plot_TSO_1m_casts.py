"""
Plots along-channel values for salinity
from Ecology CTD casts.  Uses all depths in 1 m
increments, averaged into upper and lower layers.

The goal is to see if the model does OK at reproducing
dSBAR/dx and DeltaS.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart()
import zfun
sys.path.append(os.path.abspath('../../LiveOcean/plotting'))
import pfun

# +++ load Ecology and Canadian CTD cast data +++
dir0 = Ldir['parent'] + 'ptools_data/ecology/'
sta_df = pd.read_pickle(dir0 + 'sta_df.p')

Ldir['gtagex'] = 'cas6_v3_lo8b'

TSO_all = pd.DataFrame()
TSO_dict = dict()

for year in [2017, 2018, 2019]:
    
    # input location
    in_dir = Ldir['parent'] + 'ptools_output/ecology/'
    in_fn = 'TSO_1m_casts_' + Ldir['gtagex'] + '_' + str(year) + '.p'
    
    TSO = pd.read_pickle(in_dir + in_fn)
    TSO_dict[year] = TSO
    TSO_all = pd.concat((TSO_all, TSO), sort=False)

year = 2017
this_TSO = TSO_dict[year]

out_fn = Ldir['parent'] + 'ptools_output/ecology/TSO_'+str(year)+'.png'

sn_dict = {
    0: ['SJF001','ADM002','ADM001','ADM003','PSB003','EAP001','CMB003',
         'GOR001','NSQ002','DNA001'],
    1: ['SJF001','ADM002','ADM001','HCB010','HCB003','HCB004','HCB007'],
    2: ['SJF001','ADM002','ADM001','ADM003','PSS019','SAR003']}

lab_dict = {
    0: '(a) JdF to Main Basin to South Sound',
    1: '(b) JdF to Hood Canal',
    2: '(c) JdF to Whidbey Basin'}

c_dict = {0: 'red', 1: 'olive', 2: 'orange'}

zlev_dict = {0: -20., 1: -10, 2:-10}

def make_dist(sn_list, sta_df):
    # create a distance vector
    x_list = []; y_list = []
    for sn in sn_list:
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
    return dd

# PLOTTING
fs=14
plt.rc('font', size=fs)
alpha = .5
print('\n******** YEAR = ' + str(year) + ' ********')

plt.close('all')
fig = plt.figure(figsize=(14,11))

# map
ax3 = plt.subplot2grid((3,3), (1,2), rowspan=2)
pfun.add_coast(ax3)
ax3.set_xlim(-123.5, -122)
ax3.set_ylim(47, 49)
pfun.dar(ax3)
ax3.text(.95,.05,'(d)', transform=ax3.transAxes, color='k',
    size=1.3*fs, weight='bold', ha='right')
ax3.set_xticks([-123, -122.5])
ax3.set_yticks([47,48])
ax3.set_xlabel('Longitude')
ax3.set_ylabel('Latitude')
        
for ii in [0, 1, 2]:
    
    sn_list = sn_dict[ii]
    
    zlev = zlev_dict[ii]
    
    ts_df = pd.DataFrame(index=sn_list,
            columns=['Longitude','Latitude', 'Distance (km)',
            'S top','Mod S top' ,'S bot','Mod S bot',
            'T top','Mod T top', 'T bot','Mod T bot'])

    dd = make_dist(sn_list, sta_df)
    ts_df['Distance (km)'] = dd

    for sn in sn_list:
        ts_df.loc[sn,'Longitude'] = sta_df.loc[sn,'Longitude']
        ts_df.loc[sn,'Latitude'] = sta_df.loc[sn,'Latitude']
        ts = this_TSO[this_TSO['Station']==sn]
        ts = ts[['Z','Salinity', 'Mod Salinity','Temp. (deg C)','Mod Temp. (deg C)']]
        ts = ts.dropna()
        ts_top = ts[ts['Z']>= zlev].mean()
        ts_bot = ts[ts['Z']< zlev].mean()

        ts_df.loc[sn,'S top'] = ts_top['Salinity']
        ts_df.loc[sn,'Mod S top'] = ts_top['Mod Salinity']
        ts_df.loc[sn,'S bot'] = ts_bot['Salinity']
        ts_df.loc[sn,'Mod S bot'] = ts_bot['Mod Salinity']

        ts_df.loc[sn,'T top'] = ts_top['Temp. (deg C)']
        ts_df.loc[sn,'Mod T top'] = ts_top['Mod Temp. (deg C)']
        ts_df.loc[sn,'T bot'] = ts_bot['Temp. (deg C)']
        ts_df.loc[sn,'Mod T bot'] = ts_bot['Mod Temp. (deg C)']
    
    if ii == 0:
        ax = plt.subplot2grid((3,3), (0,0), colspan=3)
        ax.text(.9,.85,'MODELED ' + str(year), weight='bold', style='italic', color='b', alpha=alpha,
            transform=ax.transAxes, ha='right', size=1.3*fs)
        ax.text(.9,.4,'OBSERVED ' + str(year), weight='bold', style='italic', color='r', alpha=alpha,
            transform=ax.transAxes, ha='right', size=1.3*fs)
    elif ii == 1:
        ax = plt.subplot2grid((3,3), (1,0), colspan=2)
    elif ii == 2:
        ax = plt.subplot2grid((3,3), (2,0), colspan=2)
        ax.set_xlabel('Distance from JdF [km]')
    ax.set_ylabel('Salinity')
            
    s0 = pd.to_numeric(ts_df['S top'].to_numpy())
    s1 = pd.to_numeric(ts_df['S bot'].to_numpy())
    
    # print stats
    print('\n'+lab_dict[ii]+', zlev = ' + str(-zlev) + ' [m]')
    BB = np.polyfit(dd,(s0+s1)/2,1)
    dsdx = BB[0] # slope
    DS = (s1-s0).mean()
    print('OBSERVED dS/dx = %0.2f [psu/100 km], DS = %0.2f [psu]' % (-100*dsdx,DS))
    ax.fill_between(dd, s0, s1, color='r', alpha=alpha)
    s0 = pd.to_numeric(ts_df['Mod S top'].to_numpy())
    s1 = pd.to_numeric(ts_df['Mod S bot'].to_numpy())
    BB = np.polyfit(dd,(s0+s1)/2,1)
    m_dsdx = BB[0] # slope
    m_DS = (s1-s0).mean()
    print('MODELED  dS/dx = %0.2f [psu/100 km], DS = %0.2f [psu]' % (-100*m_dsdx,m_DS))
    print('RATIO MOD/OBS  = %0.2f                    %0.2f' % (m_dsdx/dsdx,m_DS/DS))
    ax.fill_between(dd, s0, s1, color='b', alpha=alpha)
    ax.set_ylim(24,32)
    
    ax.text(.05,.1,lab_dict[ii], transform=ax.transAxes, color=c_dict[ii],
        size=1.3*fs, weight='bold')
        #bbox=dict(facecolor='gray', edgecolor='None', alpha=0.5))
    
    ts_df.plot(x='Longitude', y='Latitude', style='-o', color=c_dict[ii],
        linewidth=4-ii, ax=ax3, legend=False, markersize=12-(3*ii), alpha=.6)
    
plt.show()
fig.savefig(out_fn)
plt.rcdefaults()

