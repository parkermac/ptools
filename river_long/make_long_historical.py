# -*- coding: utf-8 -*-
"""
Program to gather historical records for rivers.  Focused on long
records for selected Puget Sound Rivers, relevant to the SSMSP question
of how conditions may have been different prior to 1970.
"""


#%% imports
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
import numpy as np

# SSMSP import
import os
import sys
pth = os.path.abspath('../ssmsp')
if pth not in sys.path:
    sys.path.append(pth)
import sfun
from importlib import reload
reload(sfun)

import river_class
reload(river_class)

# Load a dataframe with info for rivers to get
ri_fn = '../ssmsp/river_info.csv'

df = pd.read_csv(ri_fn, index_col='rname')

#%% set time range and select rivers

testing = False

if testing == True:
    dt0 = datetime(1900,1,1)
    dt1 = datetime(2018,12,31)
    df = df.loc[['puyallup'],:]
    save_data = False
else:
    dt0 = datetime(1900,1,1)
    dt1 = datetime(2018,12,31)
    df = df.loc[['skagit', 'snohomish', 'puyallup', 'deschutes', 'skokomish'],:]
    save_data = True
    # and create directory for output, if needed
    out_dir0 = '../../ptools_output/'
    out_dir = out_dir0 + 'river_long/'
    sfun.make_dir(out_dir0, clean=False)
    sfun.make_dir(out_dir, clean=False)

days = (dt0, dt1)

qt_dict = dict()

#%% get USGS river data

for rn in df.index:
    rs = df.loc[rn] # a series with info for this river
    riv = river_class.River(rn, rs)
    if pd.notnull(rs.usgs):
        riv.get_usgs_data(days)
        riv.print_info()
        sys.stdout.flush()
        if not riv.qt.empty:
            qt_dict[rn] = riv.qt

#%% save output

if save_data:
    for rn in qt_dict.keys():
        qt = qt_dict[rn]
        qt.to_pickle(out_dir + rn + '.p')
else:
    print('Not saving any data.')

#%% plotting

if testing:
    plt.close('all')
    fig = plt.figure(figsize=(17,9))
    ax = fig.add_subplot(111)
    qt = qt_dict['puyallup']
    qt.plot()
    plt.show()


