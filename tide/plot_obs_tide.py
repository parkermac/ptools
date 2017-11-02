"""
Code to plot year-long tide height records from
a series of NOAA and DFO sites around the Salish Sea and NE Pacific
coast.

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

# extract and save data
year  = 2013
outdir = dir0 + 'obs_data/'

M = dict()
H = dict()

for name in sn_dict.keys():
    sn = sn_dict[name]
    #fn = outdir + 'tide_' + str(sn) + '_' + str(year) + '.p'
    mfn = outdir + 'm_' + str(sn) + '_' + str(year) + '.csv'
    hfn = outdir + 'h_' + str(sn) + '_' + str(year) + '.p'
    #df = pd.read_pickle(fn)
    M[name] = Lfun.csv_to_dict(mfn)
    H[name] = pickle.load(open(hfn, 'rb'))

# plotting
plt.close('all')

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
pfun.add_coast(ax)
ax.set_xlim(-130, -122)
ax.set_ylim(42, 52)
pfun.dar(ax)

for name in sn_dict.keys():
    lon = M[name]['lon']
    lat = M[name]['lat']
    ax.plot(lon, lat, '*r')
    ax.text(lon, lat, name)
    
plt.show()
    
