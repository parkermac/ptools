# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 08:55:12 2016

@author: PM5
"""

"""
Program to gather historical TEMPERATURE records for rivers.

It gets the last 35 years for USGS rivers, but only some of them.

It gets only the last ~18 months for EC rivers.
"""

#%% imports
import matplotlib.pyplot as plt
import os
import sys

alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)

from importlib import reload
import Lfun
reload(Lfun)
Ldir = Lfun.Lstart(gridname='test')

rivp = os.path.abspath(Ldir['LO'] + 'forcing/riv1/')
if rivp not in sys.path:
    sys.path.append(rivp)

import river_T_class
reload(river_T_class)

from datetime import datetime
import pandas as pd
import numpy as np

#%% Load a dataframe with info for rivers to get

ri_fn = Ldir['parent'] + 'ptools_output/river/pnw_all_2016_07/river_info.csv'

# decide which group to get
get_usgs = False
get_ec = True

df = pd.read_csv(ri_fn, index_col='rname')

# and create directory for output, if needed
out_dir0 = Ldir['data'] + 'rivers/'
out_dir = out_dir0 + 'Data_T_historical/'
Lfun.make_dir(out_dir0, clean=False)
Lfun.make_dir(out_dir, clean=False)

#%% initialize dict

qt_dict = dict()

#%% get USGS river data

if get_usgs:

    dt0 = datetime(1980,1,1)
    dt1 = datetime(2015,12,31)
    days = (dt0, dt1)

    for rn in df.index:

        rs = df.ix[rn] # a series with info for this river
        riv = river_T_class.River(rn, rs)

        if pd.notnull(rs.usgs):# and rn in ['cedar']:

            riv.get_usgs_data(days)

            riv.print_info()
            sys.stdout.flush()

            if not riv.qt.empty:
                qt_dict[rn] = riv.qt

#%% get EC data (only has temperature for recent times)
if get_ec:

    dt0 = datetime(2015,1,1)
    dt1 = datetime(2016,6,30)
    days = (dt0, dt1)

    for rn in df.index:

        rs = df.ix[rn] # a series with info for this river

        if pd.notnull(rs.ec):# and rn in ['fraser']:

            riv = river_T_class.River(rn, rs)

            riv.get_ec_data(days)

            riv.print_info()
            sys.stdout.flush()

            if not riv.qt.empty:
                qt_dict[rn] = riv.qt

#%% save output

for rn in qt_dict.keys():

    qt = qt_dict[rn]

    qt.to_pickle(out_dir + rn + '.p')

#%% plotting

if True:

    plt.close()

    NP = len(qt_dict)

    NR = np.maximum(1, np.ceil(np.sqrt(NP)).astype(int))
    NC = np.ceil(np.sqrt(NP)).astype(int)

    fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(17,9), squeeze=False)

    cc = 0
    for rn in qt_dict.keys():
        ir = int(np.floor(cc/NC))
        ic = int(cc - NC*ir)
        ax = axes[ir, ic]
        qt_dict[rn].plot(ax=ax, title=rn, style='-k')
        ax.set_xlim(dt0, dt1)
        cc += 1

    plt.show()


