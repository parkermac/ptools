# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 09:47:51 2016

@author: PM5

This makes the file river_info.csv, and the individual river lon, lat: tracks,
for analytical (A) rivers.

"""

import os
import sys
dir0 = '/Users/PM5/Documents/'
alp = os.path.abspath(dir0 + 'LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()
import numpy as np
import pandas as pd

#%% make place for output
ri_name = 'analytical'
ri_dir0 = dir0 + 'ptools_output/river/'
ri_dir = ri_dir0 + ri_name +'/'
rit_dir = ri_dir + 'tracks/'

Lfun.make_dir(ri_dir0, clean=False)
Lfun.make_dir(ri_dir, clean=True)
Lfun.make_dir(rit_dir, clean=True)

#%% create information

names = ['creek1'] # list of river names

for rn in names:
    # save the river track information to a separate file
    if rn == 'creek1':
        N = 1000
        lon = np.linspace(-0.1, 10, N) # ndarray
        lat = 45 * np.ones(N) # ndarray
    df_tr = pd.DataFrame()
    df_tr['lon'] = lon
    df_tr['lat'] = lat
    df_tr.index.name = 'ind'
    fn_tr = rit_dir + rn + '.csv'
    df_tr.to_csv(fn_tr)

#%% initialize a DataFrame to organize the info

df = pd.DataFrame(index=names)
df['depth'] = np.nan
for rn in names:
    # save the river track information to a separate file
    if rn == 'creek1':
        df.ix[rn, 'depth'] = 5.0

#%% save to a csv file

fn_ri = ri_dir + 'river_info.csv'
df.index.name = 'rname'
df.to_csv(fn_ri)


