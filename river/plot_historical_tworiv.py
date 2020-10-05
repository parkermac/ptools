"""
Plots historical data from pickle files.

Plots two rivers at once, for the Rockfish project.

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

fs = 14
plt.rc('font', size=fs)

counter = 1
fig = plt.figure(figsize=(16,12))
for rn in ['skagit', 'fraser']:
    
    qt = pd.read_pickle(in_dir + rn + '.p')
    # NOTE: early Fraser River time index is wonky early on!
    # fixed by this sort operation
    qt = qt.sort_index()
    
    ax = fig.add_subplot(2,2,counter)

    # plot annual averages
    qa = qt.resample('A-DEC').mean()
    qa.plot(ax=ax)
    qa[[datetime(2006,12,31)]].plot(style='*r', markersize=20)
    if counter == 1:
        ax.text(.05,.15, '(a) '+rn.title()+' Annual Mean Flow', fontweight='bold', transform=ax.transAxes)
    elif counter == 2:
        ax.text(.05,.15, '(c) '+rn.title()+' Annual Mean Flow', fontweight='bold', transform=ax.transAxes)
    ax.set_ylim(0,)
    ax.grid(True)
    ax.set_xticks(pd.date_range(start='1/1/1980', end='1/1/2020', freq='5Y'))
    ax.set_xlim(datetime(1980,1,1), datetime(2020,1,1))
    ax.set_xlabel('Year')
    ax.set_ylabel('Flow $(m^{3}s^{-1})$')

    ax = fig.add_subplot(2,2,counter+2)
    
    yl = list(set(qt.index.year))
    qyd = pd.DataFrame(index=yl, columns=range(1,367))
    qyd.index.name = 'Year'
    for yr in yl:
        qy = qt[qt.index.year==yr]
        yd = qy.index.dayofyear
        qyd.loc[yr,yd] = qy.values
    # this focuses on 2006
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
    if counter==1:
        ax.set_ylim(0,5000)
        yy=3000
    elif counter==2:
        ax.set_ylim(0,14000)
        yy = 10000
    ax.fill([0,180,180,0],[0,0,yy+1000,yy+1000],'orange', alpha=.3)
    ax.fill([120,300,300,120],[0,0,yy,yy],'yellow', alpha=.3)
    ax.text(10,.8*(yy+1000),'Canary',style='italic')
    ax.text(290, .8*yy,'Yelloweye', ha='right',style='italic')
    ax.legend(loc='upper right')
    if counter == 1:
        ax.text(.05,.85, '(b) '+rn.title()+' Climatological\nFlow by Yearday', fontweight='bold', transform=ax.transAxes)
    elif counter == 2:
        ax.text(.05,.85, '(d) '+rn.title()+' Climatological\nFlow by Yearday', fontweight='bold', transform=ax.transAxes)
    ax.set_ylabel('Flow $(m^{3}s^{-1})$')
    ax.grid(True)
        
    counter += 1
    
fig.tight_layout()

plt.show()

plt.rcdefaults()
