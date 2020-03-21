"""
Plots historical data from pickle files.

"""

# imports
import matplotlib.pyplot as plt
import os, sys

sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart()

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

in_dir = Ldir['data'] + '/rivers/Data_historical/'

dt0 = datetime(2008,1,1)
dt1 = datetime(2018,12,31)

plt.close('all')

counter = 1
for rn in ['skagit', 'fraser']:#['calawah', 'hoh', 'queets', 'quinault']:
    
    if True:
        
        qt = pd.read_pickle(in_dir + rn + '.p')
        # NOTE: early Fraser River time index is wonky early on!
        # fixed by this sort operation
        qt = qt.sort_index()
        
        fig = plt.figure(figsize=(18,7))
        
        # individual plots
        # plot daily flow
        ax = fig.add_subplot(221)
        #ax = plt.subplot2grid((2,3), (0,0), colspan=2)
    
        qt.plot(ax=ax)
        ax.set_title(rn.title() + ' River')
        ax.set_ylim(0,)
        ax.grid(True)
        ax.set_xticklabels([])
        #ax.set_xlabel('')
        ax.set_xticks(pd.date_range(start='1/1/1980', end='1/1/2020', freq='5Y'))
        ax.set_xlim(datetime(1980,1,1), datetime(2020,1,1))
        ax.text(.05, .85, '(a) Daily Flow', fontweight='bold', transform=ax.transAxes)
        ax.set_ylabel('Flow $(m^{3}s^{-1})$')
    
        # plot mean flow by wateryear
        #
        # The term U.S.Geological Survey "water year" in reports that deal with surface-water
        # supply is defined as the 12-month period October 1, for any given year through
        # September 30, of the following year. The water year is designated by the calendar
        # year in which it ends and which includes 9 of the 12 months. Thus, the year ending
        # September 30, 1999 is called the "1999" water year.
        # 
        # One way to implement this would be:
        # qt.resample('A-SEP').mean()
        # which has entry:
        # 2006-09-30    2587.545205
        # and which is IDENTICAL to the result of
        # qt[datetime(2005,10,1,12):datetime(2006,9,30,12)].mean()
        
        ax = fig.add_subplot(223)
        #ax = plt.subplot2grid((2,3), (1,0), colspan=2)
        
        # NOTE, for the Skagit
        # qt.resample('AS-APR').mean()
        # Out[48]:
        # 1979-04-01    488.175415
        # 1980-04-01    460.031158
        # 1981-04-01    471.592737 **
        # 1982-04-01    509.542850
        # and
        # qt[datetime(1981,4,1,12):datetime(1982,3,31,12)].mean()
        # Out[49]: 471.592736634251
        
        if False:
            # plot by wateryear
            qa = qt.resample('A-SEP').mean()
            qa = qa[:-1] # last year is not complete
            qa.plot(ax=ax)
            qa[[datetime(2006,9,30)]].plot(style='*r', markersize=20)
            ax.text(.05,.15, '(b) Wateryear Mean Flow', fontweight='bold', transform=ax.transAxes)
            ax.set_ylim(0,)
            ax.grid(True)
            ax.set_xticks(pd.date_range(start='9/30/1980', end='9/30/2020', freq='5A-SEP'))
            ax.set_xlim(datetime(1980,9,30), datetime(2020,9,30))
            ax.set_xlabel('Wateryear')
            ax.set_ylabel('Flow $(m^{3}s^{-1})$')
        else:
            # plot annual averages
            qa = qt.resample('A-DEC').mean()
            qa.plot(ax=ax)
            qa[[datetime(2006,12,31)]].plot(style='*r', markersize=20)
            ax.text(.05,.15, '(b) Annual Mean Flow', fontweight='bold', transform=ax.transAxes)
            ax.set_ylim(0,)
            ax.grid(True)
            ax.set_xticks(pd.date_range(start='1/1/1980', end='1/1/2020', freq='5Y'))
            ax.set_xlim(datetime(1980,1,1), datetime(2020,1,1))
            ax.set_xlabel('Year')
            ax.set_ylabel('Flow $(m^{3}s^{-1})$')
    
        ax = fig.add_subplot(122)
        yl = list(set(qt.index.year))
        qyd = pd.DataFrame(index=yl, columns=range(1,367))
        qyd.index.name = 'Year'
        for yr in yl:
            qy = qt[qt.index.year==yr]
            yd = qy.index.dayofyear
            qyd.loc[yr,yd] = qy.values
        if False:
            # this bit of code focuses on long term changes
            qydm1 = qyd.loc[:2000.:].mean()
            qydm2 = qyd.loc[2000:,:].mean()
            qydm1.plot(ax=ax, label='Before 2000', style='-b')
            qydm2.plot(ax=ax, label='After 2000', style='-r')
            #
        else:
            # whereas this version focuses on 2006
            qyd_mean = qyd.mean()
            qyd_std = qyd.std()
            qyd_plus = qyd_mean + qyd_std
            qyd_minus = qyd_mean - qyd_std
            qyd_mean.plot(ax=ax, label='Daily Mean 1980-2018', style='-k', linewidth=2)
            qyd_plus.plot(ax=ax, label='+1 Std. Dev.', style='-c', linewidth=2)
            qyd_minus.plot(ax=ax, label='-1 Std. Dev.', style='-c', linewidth=2)
            qyd.loc[2006,:].plot(ax=ax, label='2006', style='-r', linewidth=3)
        ax.set_xlim(0,366)
        ax.set_xlabel('Yearday')
        ax.set_ylim(0,)
        ax.legend(loc='upper right')
        ax.text(.05,.9, '(c) Climatological\nFlow by Yearday', fontweight='bold', transform=ax.transAxes)
        ax.set_ylabel('Flow $(m^{3}s^{-1})$')
        ax.grid(True)
        
    else:
        
        if counter == 1:
            fig = plt.figure(figsize=(11,7))
        
        # all on one plot
        ax = fig.add_subplot(2,2,counter)
        qt = pd.read_pickle(in_dir + rn + '.p')
        qt.plot(ax=ax, style='-k')
        ax.set_xlim(dt0, dt1)
        ax.text(.05, .9, rn.title(), transform=ax.transAxes)
    
    counter += 1

plt.show()
