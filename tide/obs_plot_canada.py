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
indir = os.environ.get('HOME') + '/Documents/ptools_data/tide/'
name = 'tide_7795_2013.csv'
year = int(name.split('_')[-1].strip('.csv'))
fn = indir + name

df, z0 = read_tide(fn)

# and set related time limits

dt0 = datetime(year,1,1,tzinfo=pytz.timezone('UTC'))
dt1 = datetime(year+1,1,1,tzinfo=pytz.timezone('UTC'))
df = df[dt0:dt1]

dfi = read_info(fn)

# PLOTTING
plt.close('all')
figsize = (14,12)
lw = .5

fig = plt.figure(figsize=figsize)

ax = fig.add_subplot(111)
df.plot(y='z', title=('Observed Tide Height (m) ' + dfi.loc['Station_Name','Value']),
        legend=False, style='-b', ax=ax, lw=lw, grid=True, xlim=(dt0,dt1))
ax.set_xlabel('Date')
    
plt.show()

