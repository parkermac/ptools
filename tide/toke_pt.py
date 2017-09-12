"""
Code to plot observed tide time series.

Focused on Toke Pt. 2001 to explore the inverted barometer effect.

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

indir = os.environ.get('HOME') + '/Documents/ptools_data/tide/Toke_Pt_2001/'

def read_tide(in_fn):
    df = pd.read_csv(in_fn, index_col='Date Time', parse_dates = True)
    for k in df.keys():
        df = df.rename(columns={k: k.strip()})
    df = df.drop(['Sigma', 'I', 'L'], axis=1)
    df = df.rename(columns={'Water Level': 'Tide Obs'})
    # remove the mean water level
    eta0 = df['Tide Obs'].mean()
    df['Tide Obs'] = df['Tide Obs'] - eta0
    # Assumes time is UTC
    df.index.name = 'Date UTC'
    df = df.tz_localize('UTC')
    return df, eta0
    
def read_met(in_fn):
    df = pd.read_csv(in_fn, index_col='DATE TIME', parse_dates = True)
    for k in df.keys():
        df = df.rename(columns={k: k.strip()})
    df = df.drop(['RELHUM', 'VIS'], axis=1)
    # Assumes time is UTC
    df.index.name = 'Date UTC'
    df = df.tz_localize('UTC')
    return df

# READ IN OBSERVED TIDE DATA
fn = 'CO-OPS__9440910__hr.csv' # Toke Pt. 2001 observed data
obs_fn = indir + fn
obs_df, eta0 = read_tide(obs_fn)

fnm = 'CO-OPS_9440910_from_20010101_to_20011231_met.csv'
met_fn = indir + fnm
met_df = read_met(met_fn)



# make a low passed signal
eta = np.array(obs_df['Tide Obs'].tolist())
etalp = zfun.filt_godin(eta)
obs_df['Obs Low Passed'] = etalp


# merge the two
df = pd.concat([obs_df, met_df], axis=1)

# create wind stress time series
# WINDSPEED is in m/s and DIR is the compass direction
# that the wind is coming FROM [MAYBE]
#
# First: create 10m standard WSPD
P = 0.11
z_stnd = 10
z_meas = 5
df['WSPD_10'] = df['WINDSPEED'] * (z_stnd/z_meas)**P
wspd = df.WSPD_10.values
wdir = df.DIR.values
theta = 1.5*np.pi - np.pi*wdir/180.
Cd = 0.0013
rho_air = 1.22
tau = Cd * rho_air * wspd**2
taux = tau * np.cos(theta)
tauy = tau * np.sin(theta)
df['taux'] = taux
df['tauy'] = tauy
#
df['tauyF'] = df['tauy'].fillna(method='ffill')
tauy = np.array(df['tauyF'].tolist())
tauy8 = zfun.filt_AB8d(tauy)
df['NS Windstress (8-day filter)'] = tauy8


# make corrected ssh
# this is supposed to be sea level with the effect of atm pressure REMOVED
df['obs_adj'] = df['Tide Obs'] + (df['BARO'] - df['BARO'].mean())/100

# make a low passed adjusted signal
etaA = np.array(df['obs_adj'].tolist())
etaAlp = zfun.filt_godin(etaA)
df['Obs Low Passed Adj'] = etaAlp

# PLOTTING
plt.close('all')
lw = .5
fs = (14,8)

if True:
    fig = plt.figure(figsize=fs)
    ax = fig.add_subplot(311)
    df.plot(y='Tide Obs', ax=ax, lw=lw, grid=True)
    
    ax = fig.add_subplot(312)
    df.plot(y='Obs Low Passed', ax=ax, grid=True)
    df.plot(y='Obs Low Passed Adj', ax=ax, grid=True)
    
    ax = fig.add_subplot(313)
    #df.plot(y='BARO', ax=ax, grid=True)
    df.plot(y='NS Windstress (8-day filter)', ax=ax, grid=True)
    

plt.show()
