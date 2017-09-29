"""
Code to plot observed tide time series, from Canadian data.

"""

import os
import sys
import pytz
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np

from importlib import reload
import ephem_functions as efun
reload(efun)
import tractive_functions as tfun
reload(tfun)

alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import zfun

indir = os.environ.get('HOME') + '/Documents/ptools_data/tide/'

def read_tide(fn):
    df = pd.read_csv(fn, skiprows=8, header=None)
    df.columns = ['Date','z','junk']
    df = df.set_index('Date')
    df = df.drop(['junk'], axis=1)
    z0 = df['z'].mean()
    df.index = pd.to_datetime(df.index)
    df = df.tz_localize('UTC')
    return df, z0

def read_info(fn):
    df = pd.read_csv(fn, nrows=6, index_col=0, header=None)
    df.index.name = 'Item'
    df.columns = ['Value', 'junk']
    df = df.drop(['junk'], axis=1)
    return df
    
# READ IN OBSERVED TIDE DATA
fn = 'tide_7795_2013.csv'
obs_fn = indir + fn
df, z0 = read_tide(obs_fn)

dfi = read_info(obs_fn)
sname = dfi.loc['Station_Name','Value']
snum = int(dfi.loc['Station_Number','Value'])

# analysis
import utide
from matplotlib.dates import date2num
t = date2num(df.index.to_pydatetime())
z = df['z'].values

a = utide.solve(t, z, v=None,
             lat=float(dfi.loc['Latitude_Decimal_Degrees','Value']),
             nodal=False,
             trend=False,
             method='ols',
             conf_int='linear',
             Rayleigh_min=0.95)
             
plt.close('all')
fig = plt.figure()
ax = fig.add_subplot(111)
for ii in range(len(a.A)):
    # a.aux.freq has units cyles/hour
    # so for f = a.aux.frq[a.name == 'M2'][0] we get
    # 1/f = 12.420601202671868 (hours per cycle)
    ax.text(a.aux.frq[ii], a.A[ii], a.name[ii])
    
plt.show()

if False:
    dt0 = df.index[0].to_datetime()
    dt1 = df.index[-1].to_datetime()
    
    # PLOTTING
    plt.close('all')
    figsize = (14,12)
    lw = .5

    fig = plt.figure(figsize=figsize)

    ax = fig.add_subplot(111)
    df.plot(y='z', title=('Observed Tide Height (m) ' + sname),
            legend=False, style='-b', ax=ax, lw=lw, grid=True, xlim=(dt0,dt1))
    ax.set_xlabel('Date')
    
    plt.show()

