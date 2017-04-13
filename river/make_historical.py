"""
Program to gather historical records for rivers.
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

import river_class
reload(river_class)

from datetime import datetime
import pandas as pd
import numpy as np

#%% Load a dataframe with info for rivers to get

ri_fn = Ldir['parent'] + 'ptools_output/river/pnw_all_2016_07/river_info.csv'

# decide which group to get
get_usgs = False
get_ec = True
# and decide whether or not to save the data
save_data = False

df = pd.read_csv(ri_fn, index_col='rname')

# and create directory for output, if needed
out_dir0 = Ldir['data'] + 'rivers/'
out_dir = out_dir0 + 'Data_historical/'
Lfun.make_dir(out_dir0, clean=False)
Lfun.make_dir(out_dir, clean=False)

#%% set time range

dt0 = datetime(1980,1,1)
dt1 = datetime(2015,12,31)
#dt1 = datetime(1985,12,31) # debugging
days = (dt0, dt1)

qt_dict = dict()

#%% get USGS river data

if get_usgs:
    for rn in df.index:

        rs = df.ix[rn] # a series with info for this river
        riv = river_class.River(rn, rs)

        if pd.notnull(rs.usgs):

            riv.get_usgs_data(days)

            riv.print_info()
            sys.stdout.flush()

            if not riv.qt.empty:
                qt_dict[rn] = riv.qt

#%% get EC data, a year at a time
if get_ec:
    for rn in df.index:

        rs = df.ix[rn] # a series with info for this river

        Qt = pd.Series() # initialize a Series to concatenate into

        if pd.notnull(rs.ec) and rn in ['fraser']:

            #for year in range(1991, 1995): # debugging
            for year in range(dt0.year, dt1.year + 1):
                print('year = ' + str(year))

                this_days = (datetime(year,1,1), datetime(year,12,31))
                riv = river_class.River(rn, rs)

                if this_days[0] >= datetime(2015,1,1):
                    print('get 1')
                    riv.get_ec_data(this_days)
                elif this_days[0] <= datetime(2014,12,31):
                    print('get 2')
                    riv.get_ec_data_historical(year)

                riv.print_info()
                sys.stdout.flush()

                Qt = pd.concat([Qt, riv.qt])

            if not Qt.empty:
                qt_dict[rn] = Qt

#%% save output

if save_data:
    for rn in qt_dict.keys():
        qt = qt_dict[rn]
        qt.to_pickle(out_dir + rn + '.p')
else:
    print('Not saving any data.')

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


