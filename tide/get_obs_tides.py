"""
Code to automate getting year-long tide height records from
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

noaa_sn_dict = {
    'Charleston': 9432780,
    'South Beach': 9435380,
    'Garibaldi': 9437540,
    #'Cape Disappointment': 9440581,
    'Toke Point': 9440910,
    'Westport': 9441102,
    'La Push': 9442396,
    'Neah Bay': 9443090,
    'Port Angeles': 9444090,
    'Friday Harbor': 9449880,
    'Cherry Point': 9449424,
    'Port Townsend': 9444900,
    'Seattle': 9447130,
    'Tacoma': 9446484
    }
    
dfo_sn_dict = {
    'Point Atkinson': 7795,
    'Vancouver': 7735,
    'Patricia Bay': 7277,
    'Victoria Harbour': 7120,
    'Bamfield': 8545,
    'Tofino': 8615,
    'Winter Harbour': 8735,
    'Port Hardy': 8408,
    'Campbell River': 8074,
    'New Westminster': 7654
    }

# extract and save data
year  = 2013
outdir = dir0 + 'obs_data/'
Lfun.make_dir(outdir)

testing = False
if testing == True:
    noaa_sn_dict = {
        'Charleston': 9432780}
    dfo_sn_dict = {
        'Point Atkinson': 7795}

for name in noaa_sn_dict.keys():
    sn = noaa_sn_dict[name]
    fn = outdir + 'tide_' + str(sn) + '_' + str(year) + '.p'
    mfn = outdir + 'm_' + str(sn) + '_' + str(year) + '.csv'
    hfn = outdir + 'h_' + str(sn) + '_' + str(year) + '.p'
    print(name)
    df, m_dict = ofn.get_noaa_tide(sn, year)
    h = ofn.get_harmonics(df, float(m_dict['lat']))
    df.to_pickle(fn)
    Lfun.dict_to_csv(m_dict, mfn)
    pickle.dump(h, open(hfn, 'wb'))
    hh = pickle.load(open(hfn, 'rb'))

for name in dfo_sn_dict.keys():
    sn = dfo_sn_dict[name]
    fn = outdir + 'tide_' + str(sn) + '_' + str(year) + '.p'
    mfn = outdir + 'm_' + str(sn) + '_' + str(year) + '.csv'
    hfn = outdir + 'h_' + str(sn) + '_' + str(year) + '.p'
    print(name)
    df, m_dict = ofn.get_dfo_tide(sn, year)
    h = ofn.get_harmonics(df, float(m_dict['lat']))
    df.to_pickle(fn)
    Lfun.dict_to_csv(m_dict, mfn)
    pickle.dump(h, open(hfn, 'wb'))
    