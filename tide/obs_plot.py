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

# READ IN OBSERVED TIDE DATA
fn = 'CO-OPS__9447130__hr.csv' # Seattle 2016 observed data
city = 'Seattle'
obs_fn = indir + fn
obs_df, eta0 = read_tide(obs_fn)
# and set related time limits
year = 2016
dt0 = datetime(year,1,1,tzinfo=pytz.timezone('UTC'))
dt1 = datetime(year,12,31,tzinfo=pytz.timezone('UTC'))

# MAKE A SYNTHETIC TIDE
#
# General formula from Hal Mofjeld:
# h(t) = f*H*cos[ om(t-t0) + vu - G ] where
# f is the node factor and, vu = V+U is the astronomical argument
# and t0 is the start of the chosen year.
#
# load f and V+U tables from:
# http://www.pac.dfo-mpo.gc.ca/science/oceans/tidal-marees/index-eng.html
f_fn = indir + 'Foreman_node_factors.txt'
f_df = pd.read_csv(f_fn, header=0, sep='\t', index_col='Year')
vu_fn = indir + 'Foreman_astronomical_arguments.txt'
vu_df = pd.read_csv(vu_fn, header=0, sep='\t', index_col='Year')
f_ser = f_df.loc[year]
vu_ser = vu_df.loc[year]
cons_list = f_ser.keys()
#
# Load harmonic constituent values for Seattle
# Amplitudes are in meters. Phases are in degrees, referenced to GMT.
# from NOAA Tides and Currents
cons_fn = indir + 'seattle_constituents.txt'
cons_df = pd.read_csv(cons_fn, header=0, sep='\t', index_col='Name')
#
# Create the synthetic tide (inherits tzinfo from obs_df)
pred_df = pd.DataFrame(index = obs_df.index, columns=cons_list)
pred_df['tsec'] = (pred_df.index - dt0).total_seconds()
# separate constituents
for cons in cons_list:
    f = f_ser[cons]
    H = cons_df.loc[cons,'Amplitude']
    om = (np.pi/180) * cons_df.loc[cons, 'Speed']/3600
    vu = (np.pi/180) * vu_ser[cons]
    G = (np.pi/180) * cons_df.loc[cons, 'Phase']
    this_eta = f * H * np.cos(om * pred_df['tsec'] + vu - G)
    pred_df[cons] = this_eta
# create full 8-constituent prediction
pred_df['Tide Pred'] = 0
for cons in cons_list:
    pred_df['Tide Pred'] += pred_df[cons]
# create other combinations
"""
# the main diurnals:
# - can make into SunDec and MoonDec
O1                            Lunar diurnal constituent
P1                            Solar diurnal constituent
K1                            Lunar diurnal constituent (really Luni-Solar)

# the main semi-diurnals:
# - can make into Spring-Neap, with elliptical modulation
M2              Principal lunar semidiurnal constituent
S2              Principal solar semidiurnal constituent
N2        Larger lunar elliptic semidiurnal constituent

# these are small, about 5 cm:
Q1            Larger lunar elliptic diurnal constituent
K2                    Lunisolar semidiurnal constituent
"""
pred_df['SNE'] = pred_df['M2'] + pred_df['S2'] + pred_df['N2']
sun_moon_fac = ( cons_df.loc['P1','Amplitude'] /
    (cons_df.loc['P1','Amplitude'] + cons_df.loc['O1','Amplitude']) )
pred_df['SunDec'] = sun_moon_fac *pred_df['K1'] + pred_df['P1']
pred_df['MoonDec'] = (1-sun_moon_fac)*pred_df['K1'] + pred_df['O1']
pred_df['Partial Sum'] = pred_df['SNE'] + pred_df['SunDec'] + pred_df['MoonDec']

# get some orbital information
moon_orbit_df = efun.get_moon_orbit(dt0, dt1)
moon_orbit_df['Declination/20'] = moon_orbit_df['Declination (deg)']/20

# get full-new moon info
fm_df, nm_df = efun.get_full_new(dt0, dt1)
# add normalizedtractive force
fm_df = tfun.add_tf(fm_df)
nm_df = tfun.add_tf(nm_df)

# make a low passed signal
eta = np.array(obs_df['Tide Obs'].tolist())
etalp = zfun.filt_godin(eta)
obs_df['Obs Low Passed'] = etalp
obs_df['Obs-Pred'] = obs_df['Tide Obs'] - pred_df['Tide Pred']
op = np.array(obs_df['Obs-Pred'].tolist())
oplp = zfun.filt_godin(op)
obs_df['Obs-Pred Low Passed'] = oplp # appears identical to etalp

# PLOTTING
plt.close('all')
lw = .5
fs = (14,8)

if True:
    fig = plt.figure(figsize=fs)
    ax = fig.add_subplot(211)
    obs_df.plot(y='Tide Obs', ax=ax, lw=lw)
    pred_df.plot(y='Tide Pred',ax=ax, lw=lw)
    ax = fig.add_subplot(212)
    obs_df.plot(y='Obs Low Passed', ax=ax, grid=True)

if True:
    fig = plt.figure(figsize=fs)
    ax = fig.add_subplot(111)
    obs_df.plot(y='Tide Obs', ax=ax, lw=lw)
    pred_df.plot(y='SNE', ax=ax, lw=lw)
    pred_df.plot(y='SunDec', ax=ax, lw=lw)
    pred_df.plot(y='MoonDec', ax=ax, lw=lw)
    
if True:
    fig = plt.figure(figsize=fs)
    ax = fig.add_subplot(111)
    pred_df.plot(y='MoonDec', ax=ax, lw=lw)
    moon_orbit_df.plot(y='Declination/20', ax=ax)

if True:
    fig = plt.figure(figsize=fs)
    ax = fig.add_subplot(111)
    pred_df.plot(y='SNE', ax=ax, lw=lw)
    fm_df.plot(y='Tractive Force', style='or', ax=ax)
    nm_df.plot(y='Tractive Force', style='ok', ax=ax)

plt.show()
