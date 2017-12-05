"""
Code to plot observed tide time series.

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

zone='US/Pacific'
tz_local = pytz.timezone(zone)

def read_tide(in_fn):
    df = pd.read_csv(in_fn, index_col='Date Time', parse_dates = True)
    for k in df.keys():
        df = df.rename(columns={k: k.strip()})
    df = df.drop(['Sigma', 'I', 'L'], axis=1)
    df = df.rename(columns={'Water Level': 'Tide Obs'})
    # find the mean water level
    eta0 = df['Tide Obs'].mean()
    # Assumes time is UTC
    df.index.name = 'Date UTC'
    df = df.tz_localize('UTC')
    return df, eta0

# READ IN OBSERVED TIDE DATA
fn = 'CO-OPS__9447130__hr.csv' # Seattle 2016 observed data
city = 'Seattle'
obs_fn = indir + fn
obs_df, eta0 = read_tide(obs_fn)

obs_df = obs_df.tz_convert(tz_local)
obs_df.index.name = 'Date (local time)'
obs_df['Tide Obs'] = obs_df['Tide Obs'] * 3.28084

# and set related time limits
year = 2016
#tzinfo = pytz.timezone('UTC')
tzinfo = tz_local
dt0_day = datetime(year,6,10,tzinfo=tzinfo)
dt1_day = datetime(year,6,11,tzinfo=tzinfo)

dt0_month = datetime(year,6,1,tzinfo=tzinfo)
dt1_month = datetime(year,7,1,tzinfo=tzinfo)

dt0_year = datetime(year,1,1,tzinfo=tzinfo)
dt1_year = datetime(year+1,1,1,tzinfo=tzinfo)

# PLOTTING
plt.close('all')
lw0 = 0.5
lw1 = 1
lw2 = 3
fsz=18
ylim=(-5, 15)
fig = plt.figure(figsize=(14,8))

ax = fig.add_subplot(221)
obs_df.plot(y='Tide Obs',
        legend=False, style='-b', ax=ax, ylim=ylim,
        lw=lw2, grid=True, xlim=(dt0_day,dt1_day))
ax.text(.05,.05,'One Day', transform=ax.transAxes, fontweight='bold', fontsize=fsz)
ax.text(.05,.9,'Observed Tide Height (ft) ' + city,
        transform=ax.transAxes, fontsize=fsz)
ax.set_xticklabels('')
ax.set_xlabel('')

ax = fig.add_subplot(222)
obs_df.plot(y='Tide Obs',
        legend=False, style='-b', ax=ax, ylim=ylim,
        lw=lw1, grid=True, xlim=(dt0_month,dt1_month))
ax.text(.05,.05,'One Month', transform=ax.transAxes, fontweight='bold', fontsize=fsz)
ax.set_xticklabels('')
ax.set_xlabel('')

ax = fig.add_subplot(212)
obs_df.plot(y='Tide Obs',
        legend=False, style='-b', ax=ax, ylim=ylim,
        lw=lw0, grid=True, xlim=(dt0_year,dt1_year))
ax.text(.05,.05,'One Year', transform=ax.transAxes, fontweight='bold', fontsize=fsz)
ax.set_xticklabels('')
ax.set_xlabel('')

fig.set_tight_layout(True)

plt.show()

