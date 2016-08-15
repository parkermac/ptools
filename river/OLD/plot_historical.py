"""
Plots climatological data processed by river_clim.py.
"""

import matplotlib.pyplot as plt
import os; import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
from importlib import reload
import Lfun
reload(Lfun)
Ldir = Lfun.Lstart('cascadia1','base')

pth = os.path.abspath('../../LiveOcean/forcing/riv1')
if pth not in sys.path:
    sys.path.append(pth)
import river_class
reload(river_class)

import pandas as pd
import numpy as np

indir = Ldir['data'] + 'rivers/data_processed/'

if False:
    # get the list of rivers that we need for a run
    rdf = pd.read_csv(Ldir['run'] + 'rname_list.txt', header=None,
        names=['River Name'])
    rnames = rdf['River Name'].values
    rnames = rnames.tolist()
    rnames[rnames.index('duwamish')] = 'green'
    rnames[rnames.index('hammahamma')] = 'hamma'
    rnames.remove('skagit_south')
else:
    # override for testing
    rnames = [
    'columbia',
    'fraser',
    'skagit',
    'nisqually',
    'deschutes',
    'duckabush'
    ]

from datetime import datetime
dt0 = datetime(1980,1,1)
dt1 = datetime(2015,12,31)
days = (dt0, dt1)
focus_years_red = [1999] # [1998, 2002] low chinook survival
focus_years_blue = [2000] # [1999, 2004, 2005] high chinook survival

plt.close()

NR = 2; NC = 3
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(15,8), squeeze=False)

cc = 0
for rn in rnames:
    rqt = pd.read_pickle(indir + 'clim_' + rn + '.p')

    # parse the data into years, with a yearday index
    rr = pd.DataFrame(index=range(366))
    for yr in range(dt0.year, dt1.year + 1):
        qt_yr = rqt[rqt.index.year == yr]
        qt_yr.index = qt_yr.index.dayofyear
        rr[yr] = qt_yr
    rr['mean'] = rr.mean(axis=1)
    rr = rr/100.

    nr = int(np.floor(cc/NC))
    nc = int(cc - NC*nr)
    #print('cc=' + str(cc) + ', nr=' + str(nr) + ', nc=' + str(nc))
    ax = axes[nr, nc]

    for yr in range(dt0.year, dt1.year + 1):
        rr[yr].plot(ax=ax, style='-k', alpha=.2)
    rr['mean'].plot(ax=ax, style='-k', linewidth=4, alpha = .5)
    for focus_year in focus_years_red:
        rr[focus_year].plot(ax=ax, style='-r', linewidth=2)
    if len(focus_years_blue) > 0:
        for focus_year in focus_years_blue:
            rr[focus_year].plot(ax=ax, style='-b', linewidth=2)
    if nc == 0:
        ax.set_ylabel('Flow (100 $m^{3}s^{-1}$)')
    if nr == NR-1:
        ax.set_xlabel('Yearday')
    else:
        ax.xaxis.set_ticklabels([])
    ax.set_xlim(rr.index[0], rr.index[-1])
    isgood = rr.ix[160][:].notnull()
    first_good_year = rr.columns[isgood][0]
    title_text =  rn.title() + ' ' + str(first_good_year) + '-' + str(dt1.year)
    ax.text(0.05, 0.9, title_text, transform=ax.transAxes)
    if nr == 0 and nc == 0:
        ax.text(0.95, 0.9, 'Mean', horizontalalignment='right',
            color='k', alpha=0.5, transform=ax.transAxes)
        ax.text(0.95, 0.85, str(focus_years_red), horizontalalignment='right',
            color='r', transform=ax.transAxes)
        if len(focus_years_blue) > 0:
            ax.text(0.95, 0.8, str(focus_years_blue), horizontalalignment='right',
                color='b', transform=ax.transAxes)
    ax.grid()

    cc += 1

plt.show()