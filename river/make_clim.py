# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 13:50:31 2016

@author: PM5

Program to create climatological records for rivers.
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

import pandas as pd

#%% Load a dataframe with info for rivers to get

ri_fn = Ldir['parent'] + 'ptools_output/river/pnw_all_2016_07/river_info.csv'
df = pd.read_csv(ri_fn, index_col='rname')

# create directory for output, if needed
dir0 = Ldir['data'] + 'rivers/'
in_dir = dir0 + 'Data_historical/'
out_dir = dir0 + 'Data_clim/'
Lfun.make_dir(dir0, clean=False)
Lfun.make_dir(out_dir, clean=False)

#%% create the climatologies

for rn in df.index:

    qt = pd.read_pickle(in_dir + rn + '.p')

    clm_d = dict()

    for day in range(1,367):
        clm_d[day] = qt[qt.index.dayofyear == day].mean()

    qtc = pd.Series(clm_d)

    # plot if desired
    if rn in ['skagit']:
        fig = plt.figure(figsize=(15,10))
        ax = fig.add_subplot(111)
        year0 = qt.index[0].year
        year1 = qt.index[-1].year
        for year in range(year0, year1 + 1):
            this_qt = qt[qt.index.year == year]
            ax.plot(this_qt.index.dayofyear, this_qt.values, alpha=.3)
        ax.plot(qtc.index, qtc.values, color='r', linewidth=4)
        ax.set_title(rn.title())
        ax.set_xlim(0,366)
        ax.set_xlabel('Yearday')

    # save climatology to a csv file
    qtc.to_csv(out_dir + rn + '.csv')