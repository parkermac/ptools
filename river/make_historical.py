"""
Program to gather historical records for rivers.
"""

#%% imports
import matplotlib.pyplot as plt
import os
import sys
from importlib import reload

alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()
import zfun
reload(zfun)

rivp = os.path.abspath(Ldir['LO'] + 'forcing/riv2/')
if rivp not in sys.path:
    sys.path.append(rivp)
import river_class
reload(river_class)

from datetime import datetime, timedelta
import pandas as pd
import numpy as np

#%% Load a dataframe with info for rivers to get
ri_fn = Ldir['parent'] + 'ptools_output/river/pnw_all_2016_07/river_info.csv'

# decide which group to get
get_usgs = True
get_ec = True
# and decide whether or not to save the data
save_data = True

df = pd.read_csv(ri_fn, index_col='rname')


#%% set time range

testing = False # custom settings

if testing == True:
    dt0 = datetime(2013,1,1)
    dt1 = datetime(2019,2,14)
    #df = df.loc[['skokomish', 'nf_skokomish', 'sf_skokomish', 'fraser']]
    #df = df.loc[['columbia', 'naselle', 'willapa']]
    df = df.loc[['fraser']]
    save_data = True
    # and create directory for output, if needed
    out_dir00 = Ldir['parent'] + 'ptools_output/'
    out_dir0 = out_dir00 + 'river/'
    out_dir = out_dir0 + 'ec_testing_2019.02.17/'
    Lfun.make_dir(out_dir00, clean=False)
    Lfun.make_dir(out_dir0, clean=False)
    Lfun.make_dir(out_dir, clean=False)

else:
    dt0 = datetime(1980,1,1)
    dt1 = datetime(2018,12,31)
    # and create directory for output, if needed
    out_dir0 = Ldir['data'] + 'rivers/'
    out_dir = out_dir0 + 'Data_historical/'
    Lfun.make_dir(out_dir0, clean=False)
    Lfun.make_dir(out_dir, clean=False)
days = (dt0, dt1)

qt_dict = dict()

#%% get USGS river data

if get_usgs:
    for rn in df.index:
        rs = df.loc[rn] # a series with info for this river
        riv = river_class.River(rn, rs)
        if pd.notnull(rs.usgs):
            riv.get_usgs_data(days)
            riv.print_info()
            sys.stdout.flush()
            if not riv.qt.empty:
                qt_dict[rn] = riv.qt

#%% get EC data, a year at a time
roms_fn = 'cas4_v2_2017.01.01_2018.12.31.p'
roms_qt = pd.read_pickle(Ldir['LOo'] + 'river/' + roms_fn)
if get_ec:
    for rn in df.index:
        rs = df.loc[rn] # a series with info for this river
        Qt = pd.Series() # initialize a Series to concatenate into
        if pd.notnull(rs.ec) and rn in ['fraser']:
            for year in range(dt0.year, dt1.year + 1):
                print('year = ' + str(year))
                this_days = (datetime(year,1,1), datetime(year,12,31))
                riv = river_class.River(rn, rs)
                if year >= 2019:
                    print(' - getting current EC data')
                    riv.get_ec_data(this_days)
                    riv.print_info()
                    sys.stdout.flush()
                    this_qt = riv.qt
                elif year in [2017, 2018]:
                    this_qt = roms_qt.loc[this_days[0]:this_days[1], rn]
                    this_qt.index = this_qt.index + timedelta(days=0.5)
                    print(' - getting historical EC data from ' + roms_fn)
                elif year <= 2016:
                    print(' - getting historical EC data')
                    riv.get_ec_data_historical(year)
                    riv.print_info()
                    sys.stdout.flush()
                    this_qt = riv.qt
                Qt = pd.concat([Qt, this_qt])
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

    plt.close('all')

    NP = len(qt_dict)
    NR, NC = zfun.get_rc(NP)
    fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(17,9), squeeze=False)
    ii = 0
    for rn in qt_dict.keys():
        ir, ic = zfun.get_irc(ii, NC)
        ax = axes[ir, ic]
        this_ser = qt_dict[rn]
        this_ser.plot(ax=ax, style='-k')
        ax.set_xlim(dt0, dt1)
        ax.text(.05, .9, rn, transform=ax.transAxes)
        ii += 1

    plt.show()


