"""
Plot low-passed SSH fields from model and observations for
selected stations.

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
year  = 2013

sn_list = ['Neah Bay', 'Vancouver', 'Seattle']

obs_dir = dir0 + 'obs_data/'
Tobs = dict()
Mobs = dict()
Hobs = dict()
for name in sn_list:
    # observed tide
    sn = sn_dict[name]
    fn = obs_dir + 'tide_' + str(sn) + '_' + str(year) + '.p'
    mfn = obs_dir + 'm_' + str(sn) + '_' + str(year) + '.csv'
    hfn = obs_dir + 'h_' + str(sn) + '_' + str(year) + '.p'
    Tobs[name] = pd.read_pickle(fn)
    Mobs[name] = Lfun.csv_to_dict(mfn)
    Hobs[name] = pickle.load(open(hfn, 'rb'))
    
mod_dir = dir0 + 'mod_data/'
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
ax.set_ylim(44, 50)
pfun.dar(ax)
for name in sn_list:
    lon = float(Mobs[name]['lon'])
    lat = float(Mobs[name]['lat'])
    ax.plot(lon, lat, '*r')
    ax.text(lon+.03, lat, name, verticalalignment='center')
    
def add_lp(df):
    # add low-passed signal
    eta = np.array(df['eta'].tolist())
    eta0 = np.mean(eta)
    eta_lp = zfun.filt_godin(eta - eta0)
    df['eta_lp'] = eta_lp
    return df, eta0
    
# data panels

# note that we are not correcting for atmospheric pressure

NS = len(sn_list)
for ii in range(NS):
    name = sn_list[ii]
    ax = plt.subplot2grid((NS,3), (ii,0), colspan=2)
    Tm = Tmod[name]
    To = Tobs[name]
    Tm, eta0m = add_lp(Tm)
    To, eta0o = add_lp(To)
    Tm.plot(y='eta_lp', style='-b', ax=ax, legend=False)
    To.plot(y='eta_lp', style='-r', ax=ax, legend=False, ylim=(-.4, .4), grid=True)
    ax.text(.05, .9,
        name,
        transform=ax.transAxes, weight='bold')
    if ii < NS-1:
        ax.set_xticklabels('')
        ax.set_xlabel('')
        
plt.show()
    
