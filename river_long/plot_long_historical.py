# -*- coding: utf-8 -*-
"""
Program to plot historical records for rivers.  Focused on long
records for selected Puget Sound Rivers, relevant to the SSMSP question
of how conditions may have been different prior to 1970.
"""


# imports
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime, timedelta

out_dir0 = '../../ptools_output/'

# Load a dataframe with info for selected rivers
ri_fn = '../ssmsp/river_info.csv'
df = pd.read_csv(ri_fn, index_col='rname')
df = df.loc[['skagit', 'snohomish', 'puyallup', 'deschutes', 'skokomish'],:]
#df = df.loc[['skagit'],:]

# where the extracted data is
out_dir = out_dir0 + 'river_long/'

plt.close('all')

for rn in df.index:
    fig = plt.figure(figsize=(16,8))

    qt = pd.read_pickle(out_dir + rn + '.p')
    
    # drop bad data from later years of Skokomish - they stopped reporting high flow events
    if rn == 'skokomish':
        qt = qt[qt.index[0]:pd.datetime(2008,12,31)]

    # plot daily flow
    ax = fig.add_subplot(221)
    #ax = plt.subplot2grid((2,3), (0,0), colspan=2)
    
    qt.plot(ax=ax)
    ax.set_title(rn.title() + ' River')
    ax.set_ylim(0,)
    ax.grid(True)
    ax.set_xticklabels([])
    #ax.set_xlabel('')
    ax.set_xticks(pd.date_range(start='1/1/1900', end='1/1/2020', freq='10Y'))
    ax.set_xlim(qt.index[0], pd.datetime(2020,1,1))
    #ax.set_xlim(pd.datetime(1900,1,1), pd.datetime(2020,1,1))
    ax.text(.05, .85, '(a) Daily Flow', fontweight='bold', transform=ax.transAxes)
    ax.set_ylabel('Flow $(m^{3}s^{-1})$')
    
    # plot mean flow by wateryear
    ax = fig.add_subplot(223)
    #ax = plt.subplot2grid((2,3), (1,0), colspan=2)
    
    wyt = qt.index + timedelta(days=92)
    qwy = pd.Series(index=wyt, data=qt.values)
    qa = qwy.resample('Y').mean()
    qa = qa[:-1] # drop last wateryear (2019) because it is not over
    
    qa.plot(ax=ax)
    ax.set_ylim(0,)
    ax.grid(True)
    ax.set_xticks(pd.date_range(start='1/1/1900', end='1/1/2020', freq='10Y'))    
    ax.set_xlim(qt.index[0], pd.datetime(2020,1,1))
    #ax.set_xlim(pd.datetime(1900,1,1), pd.datetime(2020,1,1))
    ax.set_xlabel('Year')
    ax.text(.05,.15, '(b) Annual Mean Flow by Wateryear', fontweight='bold', transform=ax.transAxes)
    ax.set_ylabel('Flow $(m^{3}s^{-1})$')
    
    ax = fig.add_subplot(122)
    yl = list(set(qt.index.year))
    qyd = pd.DataFrame(index=yl, columns=range(1,367))
    qyd.index.name = 'Year'
    for yr in yl:
        qy = qt[qt.index.year==yr]
        yd = qy.index.dayofyear
        qyd.loc[yr,yd] = qy.values
    qydm1 = qyd.loc[:1970.:].mean()
    qydm2 = qyd.loc[1970:,:].mean()
    qydm1.plot(ax=ax, label='Before 1970', style='-b')
    qydm2.plot(ax=ax, label='After 1970', style='-r')
    ax.set_xlim(0,366)
    ax.set_xlabel('Yearday')
    ax.set_ylim(0,)
    ax.legend(loc='lower left')
    ax.text(.05,.9, '(c) Climatological Flow by Yearday', fontweight='bold', transform=ax.transAxes)
    ax.set_ylabel('Flow $(m^{3}s^{-1})$')
    ax.grid(True)
    
    plt.savefig(out_dir + 'Historical_' + rn + '.png')
    
plt.show()
    
