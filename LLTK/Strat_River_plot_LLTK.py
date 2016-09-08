"""
Compare stratification and river flow.

May be obsolete.  2016.09.08 PM
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

#%% STRATIFICATION
# These are saved as Series of Sigma (top or bottom layer)

# where are the data csv files
dir0 = '/Users/PM5/Documents/'
indir0 = dir0 + 'ptools_output/env_data/'
indir = indir0 + 'EcologyCTD/'

# choose which stations to include by commenting out parts of this list
# names ending in _0 are multi-year 1990-2014
sta_to_plot = [
    'BLL009_0', # Bellingham Bay
    'BUD005_0', # Budd Inlet
    'DNA001_0', # Dana Passage
                # 2001 deep T bad in July
                # 2002 shallow S low value Feb, shallow T spike Oct
                # 2006 lots of little T, s spikes all depths, many months
    'HCB004_0', # Hood Canal in Lynch Cove
    'PSB003_0', # Main Basin off Seattle
                # 1994 salinity spike 4 m Feb/Mar
                # 2007 many T, s spikes, especially deep
    'HCB003_0',   # Hood Canal, middle of main channel (Hama Hama)
    'SAR003_0',   # middle of Whidbey Basin
    'GOR001_0',
    'PSS019_0',
    'PSB003_0',
    'CRR001_0',
    ]

ser_top_dict = dict()
ser_bot_dict = dict()
ser_strat_dict = dict()
for sta in sta_to_plot:
    infile = indir + 'ser_top_' + sta + '.p'
    ser_top_dict[sta] = pd.read_pickle(infile)
    infile = indir + 'ser_bot_' + sta + '.p'
    ser_bot_dict[sta] = pd.read_pickle(infile)
    ser_strat_dict[sta] = ser_bot_dict[sta] - ser_top_dict[sta]

def month_year_bins(ser):
    s = ser
    # parse the data into years, with a month index
    cols = range(1990, 2016)
    s_clim = pd.DataFrame(index=range(1,13), columns=cols)
    for yr in range(1990, 2015 + 1):
        sy = s[s.index.year == yr]
        for m in range(1,13):
            sym = sy[sy.index.month == m].values
            if len(sym) == 0:
                symm = np.nan
            else:
                symm = sym.mean()
            s_clim.ix[m][yr] = symm
    s_clim['mean'] = s_clim.mean(axis=1)
    return s_clim

clim_top_dict = dict()
clim_bot_dict = dict()
clim_strat_dict = dict()

for sta in sta_to_plot:
    clim_top_dict[sta] = month_year_bins(ser_top_dict[sta])
    clim_bot_dict[sta] = month_year_bins(ser_bot_dict[sta])
    clim_strat_dict[sta] = month_year_bins(ser_strat_dict[sta])

#%% RIVERS
# these are saved in DataFrames with year columns (and mean)

pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart('cascadia1','base')

# get the list of rivers that we need for a run
rdf = pd.read_csv(Ldir['grid'] + 'rname_list.txt', header=None,
    names=['River Name'])
rnames = rdf['River Name'].values
rnames = rnames.tolist()
rnames[rnames.index('duwamish')] = 'green'
rnames[rnames.index('hammahamma')] = 'hamma'
rnames.remove('skagit_south')

indir = Ldir['data'] + 'rivers/data_processed/'

riv_dict = dict()
clim_riv_dict = dict()
for rn in rnames:
    rqt = pd.read_pickle(indir + 'clim_' + rn + '.p')
    riv_dict[rn] = rqt
    clim_riv_dict[rn] = month_year_bins(rqt)

sta_to_plot = [
    #'BLL009_0',
    #'BUD005_0',
    #'DNA001_0',
    'HCB004_0',
    #'PSB003_0',
    #'HCB003_0',
    #'SAR003_0',
    'GOR001_0',
    'PSS019_0',
    'PSB003_0',
    'CRR001_0',
    ]
riv_for_sta = [
    #'fraser',
    #'deschutes',
    #'deschutes',
    'skokomish',
    #'puyallup',
    #'dosewallips',
    #'skagit',
    'nisqually',
    'snohomish',
    'green',
    'nisqually',
    ]
s2r = dict(zip(sta_to_plot,riv_for_sta))

#%% make aligned time series for all river-station pairs
aligned_strat_ser = dict()
aligned_riv_ser = dict()
for sta in sta_to_plot:

    # stratification
    strat_dict = ser_strat_dict[sta]
    strat_ser = pd.Series(strat_dict['Sigma'], index=strat_dict.index)
    new_ind = []
    for ssi in strat_ser.index:
        new_ind.append(datetime(ssi.year,ssi.month,1))
    strat_ser.index = new_ind
    strat_ser = strat_ser.dropna()
    # this drops entries with duplicated indices
    strat_ser = strat_ser.groupby(strat_ser.index).first()

    # rivers
    rn = s2r[sta]
    riv_ser = riv_dict[rn]
    riv_ser.index.name = 'date'
    # this gives Qr values averaged over the past month,
    # with the index at the first of the month
    # (e.g. index = February 1 for average of January)
    riv_ser = riv_ser.resample('M', how='mean')
    riv_ser.index = riv_ser.index + timedelta(1)
    riv_ser = riv_ser.dropna()
    riv_ser = riv_ser.groupby(riv_ser.index).first()

    # combining indices
    joint_ind = riv_ser.index.intersection(strat_ser.index)
    #stix = strat_ser[joint_ind]
    #riix = riv_ser[joint_ind]
    stix = pd.Series()
    riix = pd.Series()
    for item in joint_ind:
        stix[item] = strat_ser[item]
        riix[item] = riv_ser[item]
    print('River: %s Station: %s ji(%d) stix(%d) riix(%d)' %
        (rn, sta, len(joint_ind), len(stix), len(riix)))
    aligned_strat_ser[sta] = stix
    aligned_riv_ser[sta] = riix

#%% plotting
plt.close('all')

NR = 3
NC = len(sta_to_plot)
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(15,10), squeeze=False)
counter = 0
for sta in sta_to_plot:
    rn = s2r[sta]
    ax0 = axes[0][counter]
    ax1 = axes[1][counter]
    ax2 = axes[2][counter]
    s = clim_strat_dict[sta]['mean']
    r = clim_riv_dict[rn]['mean']
    sm = s.values
    rm = r.values
    ax0.plot((s.index.values - .5), sm, '-b', linewidth=3)
    ax0.set_title(sta)
    if counter == 0:
        ax0.set_ylabel('Mean Strat $(kg m^{-3})$')
        ax1.set_ylabel('Mean Flow $(100 m^{3} s^{-1})$')
        ax2.set_ylabel('Strat $(kg m^{-3})$')
    ax1.set_title(rn.title())
    ax1.set_xlabel('Month')
    ax1.plot((r.index.values - .5), rm/100, '-b', linewidth=3)
    ax0.set_xlim(0,12)
    ax1.set_xlim(0,12)

    sv = aligned_strat_ser[sta].values
    rv = aligned_riv_ser[sta].values

    lag_list = range(-2,5)
    cc_list = []
    #print('\n** River = ' + rn + '\n** Station = ' + sta)
    for lag in lag_list:
        cc = np.corrcoef(np.roll(rv,lag), sv)[0,1]
        cc_list.append(cc)
        #print(' lag = %s months: cc = %0.2f' % (lag, cc))

    n_bestlag = np.argmax(cc_list)
    bestlag = lag_list[n_bestlag]
    ax2.plot(np.roll(rv,bestlag)/100, sv, '*r')
    ax2.text(.05,.8,'lag = %d mo\n cc = %0.2f' % (bestlag, cc_list[n_bestlag]),
        transform=ax2.transAxes)
    ax2.set_xlabel('Flow $(100 m^{3} s^{-1})$')

    counter += 1

#NR = 1
#NC = len(sta_to_plot)
#fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(15,8), squeeze=False)
#counter = 0
#for sta in sta_to_plot:
#    rn = s2r[sta]
#    ax0 = axes[0][counter]
#    astr = aligned_strat_ser[sta]
#    ariv = aligned_riv_ser[sta]
#    ax0.plot(ariv.values, astr.values, '*r')
#    ax0.set_title(sta)
#    cc = np.corrcoef(ariv.values, astr.values)[0,1]
#    ax0.text(.05,.8,'cc = %0.2f' % (cc), transform=ax0.transAxes)
#
#    counter += 1

plt.show()



