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

# conversion factors
m2f = 3.28
f2m = 1/3.28

indir = os.environ.get('HOME') + '/Documents/ptools_data/tide/Toke_Pt_2001/'

def read_tide(in_fn):
    df = pd.read_csv(in_fn, index_col='Date Time', parse_dates = True)
    for k in df.keys():
        df = df.rename(columns={k: k.strip()})
    df = df.drop(['Sigma', 'I', 'L'], axis=1)
    df = df.rename(columns={'Water Level': 'z'})
    # find the mean water level
    z0 = df['z'].mean()
    # Assumes time is UTC
    df.index.name = 'Date'
    df = df.tz_localize('UTC')
    return df, z0
    
def read_met(in_fn):
    df = pd.read_csv(in_fn, index_col='DATE TIME', parse_dates = True)
    for k in df.keys():
        df = df.rename(columns={k: k.strip()})
    df = df.drop(['RELHUM', 'VIS'], axis=1)
    # Assumes time is UTC
    df.index.name = 'Date'
    df = df.tz_localize('UTC')
    return df

# Load tide data
fn = 'CO-OPS__9440910__hr.csv' # Toke Pt. 2001 observed data
obs_fn = indir + fn
obs_df, z0 = read_tide(obs_fn)

# Load met data
fnm = 'CO-OPS_9440910_from_20010101_to_20011231_met.csv'
met_fn = indir + fnm
met_df = read_met(met_fn)

# merge the two
df = pd.concat([obs_df, met_df], axis=1)

# create wind stress time series
# WINDSPEED is in m/s and DIR is the compass direction
# that the wind is coming FROM [MAYBE]
#
# First: create 10m standard WSPD
P = 0.11
z_stnd = 10
z_meas = 22.6 * f2m
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
# fill in gaps
df['tauy'] = df['tauy'].fillna(method='ffill')
#
# filter the windstress
tauy = np.array(df['tauy'].tolist())
tauy8 = zfun.filt_AB8d(tauy)
df['tauy8'] = tauy8

# make adjusted ssh
# this is supposed to be sea level with the effect of atm pressure REMOVED
# NOTE: pressure is in millibars
df['BARO'] = df['BARO'].fillna(method='ffill')# fill gaps
df['p'] = -(df['BARO'] - df['BARO'].mean())/100
df['za'] = df['z'] + (df['BARO'] - df['BARO'].mean())/100

# make low-passed signal
z = np.array(df['z'].tolist())
zlp = zfun.filt_godin(z)
df['zlp'] = zlp

# make low-passed adjusted signal
za = np.array(df['za'].tolist())
zalp = zfun.filt_godin(za)
df['zalp'] = zalp

# make low-passed atm pressure contribution
p = np.array(df['p'].tolist())
plp = zfun.filt_godin(p)
df['plp'] = plp

# convert to feet
df['z'] = df['z'] * m2f
df['zlp'] = df['zlp'] * m2f
df['za'] = df['za'] * m2f
df['zalp'] = df['zalp'] * m2f
df['plp'] = df['plp'] * m2f

# make separate positive and negative windstress
tauy8p = np.array(df['tauy8'].tolist())
tauy8p[np.isnan(tauy8p)] = 0
tauy8p[tauy8p < 0] = 0
df['tauy8p'] = tauy8p

tauy8m = np.array(df['tauy8'].tolist())
tauy8m[np.isnan(tauy8m)] = 0
tauy8m[tauy8m > 0] = 0
df['tauy8m'] = tauy8m

# PLOTTING
plt.close('all')
fs = (14,8)

fsz = 14

if True:
    fig = plt.figure(figsize=fs)
    
    ax = fig.add_subplot(311)
    df.plot(y='z',
        lw=0.5, color='k',
         ax=ax, grid=True, legend=False)
    ax.set_xticklabels('')
    ax.set_xlabel('')
    ax.text(.05,.85,'Toke Pt. Tide Height (feet)', color='k',
        transform=ax.transAxes, fontweight='bold', fontsize=fsz)
    
    
    ax1 = fig.add_subplot(312)
            
    ax1.fill_between(df.index, df['tauy8p'].values, color='r')
    ax1.fill_between(df.index, df['tauy8m'].values, color='b')
    
    ax1.grid()
    ax1.set_xlim((df.index[0], df.index[-1]))
    
    ax1.set_xticklabels('')
    ax1.set_xlabel('')
    
    ax1.text(.05,.05,'North-South Wind UPWELLING', color='b',
        transform=ax1.transAxes, fontweight='bold', fontsize=fsz)
    ax1.text(.05,.85,'North-South Wind DOWNWELLING', color='r',
        transform=ax1.transAxes, fontweight='bold', fontsize=fsz)
    
    ax2 = fig.add_subplot(313)
    
    df.plot(y='zlp',
        lw=3, color='k',
        ax=ax2, grid=True, legend=False)
        
    # df.plot(y='p',
    #     lw=2, color='g',
    #     ax=ax2, grid=True, legend=False)
    
    df.plot(y='zalp',
        lw=2, color='g',
        ax=ax2, grid=True, legend=False)
        
    ax2.text(.05,.85,'Daily Mean Tide Height (ft)', color='k',
        transform=ax2.transAxes, fontweight='bold', fontsize=fsz)
    ax2.text(.05,.7,'Daily Mean Tide Height (ft) Without Atmospheric Pressure', color='g',
        transform=ax2.transAxes, fontweight='bold', fontsize=fsz)

        
    ax2.set_xlabel('Date')
    

plt.show()
