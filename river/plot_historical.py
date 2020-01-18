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
for rn in ['calawah', 'hoh', 'queets', 'quinault']:
    
    if True:
        
        qt = pd.read_pickle(in_dir + rn + '.p')
        
        fig = plt.figure(figsize=(11,7))
        
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
        ax.set_xlim(pd.datetime(1980,1,1), pd.datetime(2020,1,1))
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
        ax.set_xticks(pd.date_range(start='1/1/1980', end='1/1/2020', freq='5Y'))
        ax.set_xlim(pd.datetime(1980,1,1), pd.datetime(2020,1,1))
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
        qydm1 = qyd.loc[:2000.:].mean()
        qydm2 = qyd.loc[2000:,:].mean()
        qydm1.plot(ax=ax, label='Before 2000', style='-b')
        qydm2.plot(ax=ax, label='After 2000', style='-r')
        ax.set_xlim(0,366)
        ax.set_xlabel('Yearday')
        ax.set_ylim(0,)
        ax.legend(loc='lower left')
        ax.text(.05,.9, '(c) Climatological Flow by Yearday', fontweight='bold', transform=ax.transAxes)
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
