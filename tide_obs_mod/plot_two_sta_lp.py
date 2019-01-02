"""
Plot low-passed SSH fields from model for a pair of stations.

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import zfun
import zrfun

pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np
import pickle

from importlib import reload
import obsfun as ofn
reload(ofn)

home = os.environ.get('HOME')
dir00 = home + '/Documents/'
dir0 = dir00 + 'ptools_output/tide/'

noaa_sn_dict, dfo_sn_dict, sn_dict = ofn.get_sn_dicts()

# load data
year  = 2017

sn_list = ['Neah Bay', 'Campbell River']

    
mod_dir = dir0 + 'mod_data/cas4_v2_lo6biom/'
Tmod = dict()
Mmod = dict()
Hmod = dict()
for name in sn_list:
    # observed tide
    sn = sn_dict[name]
    fn = mod_dir + 'tide_' + str(sn) + '_' + str(year) + '.p'
    mfn = mod_dir + 'm_' + str(sn) + '_' + str(year) + '.csv'
    hfn = mod_dir + 'h_' + str(sn) + '_' + str(year) + '.p'
    Tmod[name] = pd.read_pickle(fn)
    Mmod[name] = Lfun.csv_to_dict(mfn)
    Hmod[name] = pickle.load(open(hfn, 'rb'))

# plotting
plt.close('all')
fig = plt.figure(figsize=(14,8))

# station map
ax = plt.subplot2grid((1,3), (0,2), colspan=1)
pfun.add_coast(ax)
ax.set_xlim(-126, -122)
ax.set_ylim(44, 51)
pfun.dar(ax)
for name in sn_list:
    lon = float(Mmod[name]['lon'])
    lat = float(Mmod[name]['lat'])
    ax.plot(lon, lat, '*r')
    ax.text(lon+.03, lat, name, verticalalignment='center')
    
def add_lp(df):
    # add low-passed signal
    eta = np.array(df['eta'].tolist())
    eta_lp = zfun.filt_godin(eta)
    df['eta_lp'] = eta_lp
    return df
    
df = pd.DataFrame(index=Tmod[sn_list[0]].index, columns = sn_list + ['dz'])

for name in sn_list:
    Tm = Tmod[name]
    Tm= add_lp(Tm)
    df[name] = Tm['eta_lp']
df['dz'] = df[sn_list[1]] - df[sn_list[0]]

ax = plt.subplot2grid((2,3), (0,0), colspan=2)
df.plot(y=sn_list, ax=ax, legend=True, grid=True)
ax.set_xticklabels('')
ax.set_xlabel('')
ax.set_xlim(datetime(2017,1,1), datetime(2018,1,1))

ax = plt.subplot2grid((2,3), (1,0), colspan=2)
df.plot(y='dz', ax=ax, legend=True, grid=True)
ax.set_xlim(datetime(2017,1,1), datetime(2018,1,1))

plt.show()
    
