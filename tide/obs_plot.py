"""
Code to plot observed tide time series.

"""

import os
# import sys
# pth = os.path.abspath('../../LiveOcean/alpha')
# if pth not in sys.path:
#     sys.path.append(pth)
# import Lfun

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np

# from importlib import reload
# import transit
# reload(transit)

def read_tide(in_fn, city):
    df = pd.read_csv(in_fn, index_col='Date Time', parse_dates = True)
    for k in df.keys():
        df = df.rename(columns={k: k.strip()})
    df = df.drop(['Sigma', 'I', 'L'], axis=1)
    df = df.rename(columns={'Water Level': 'Tide ' + city})
    df.index.name = 'Date UTC'
    return df

home = os.environ.get('HOME')
dir00 = home + '/Documents/'

# read in tide data (time is UTC)
indir = dir00 + 'ptools_data/tide/'

# Seattle 2016 observed data
fn = 'CO-OPS__9447130__hr.csv'
city = 'Seattle'
df1 = read_tide(indir+fn, city)

# fn = 'CO-OPS__9441102__hr.csv'
# city = 'Westport'
# df2 = read_tide(indir+fn, city)

# make a synthetic tide
#
# General formula from Hal Mofjeld
# h(t) = f*H*cos[ om(t-t0) + vu - G ]
#
# Nodal corrections, 2016 values, from:
# http://www.pac.dfo-mpo.gc.ca/science/oceans/tidal-marees/index-eng.html
t0 = datetime(2016,1,1,0,0)
cons_list = ['Q1','O1','P1','K1','N2','M2','S2','K2']
f_list = [0.821,0.809,1.011,0.886,1.037,1.037,0.998,0.752]
vu_list = [43.05,202.05,349.847,9.242,50.69,210.745,0.001,198.451]
f_dict = dict(zip(cons_list, f_list))
vu_dict = dict(zip(cons_list, vu_list))
#
# Constituent values for Seattle
# Amplitudes are in meters. Phases are in degrees, referenced to GMT.
# from NOAA Tides and Currents
c_fn = indir + 'seattle_constituents.txt'
df_cons = df = pd.read_csv(c_fn, header=0, sep='\t', index_col='Name')
#
# Create the synthetic tide
df3 = pd.DataFrame(index = df1.index, columns=cons_list)
df3['tsec'] = (df3.index - t0).total_seconds()
# separate constituents
for cons in cons_list:
    f = f_dict[cons]
    H = df_cons.loc[cons,'Amplitude']
    om = (np.pi/180) * df_cons.loc[cons, 'Speed']/3600
    vu = (np.pi/180) * vu_dict[cons]
    G = (np.pi/180) * df_cons.loc[cons, 'Phase']
    this_eta = f * H * np.cos(om * df3['tsec'] + vu - G)
    df3[cons] = this_eta
# full 8-constituent prediction
df3['Tide Prediction'] = 0
for cons in cons_list:
    df3['Tide Prediction'] += df3[cons]
eta0 = float(df1['Tide Seattle'].mean())
df3['Tide Prediction'] += eta0
# other combinations
"""
# the main diurnals:
# - can make into SunDec and MoonDec
O1                            Lunar diurnal constituent
P1                            Solar diurnal constituent
K1                            Lunar diurnal constituent

# the main semi-diurnals:
# - can make into Spring-Neap, with elliptical modulation
M2              Principal lunar semidiurnal constituent
S2              Principal solar semidiurnal constituent
N2        Larger lunar elliptic semidiurnal constituent

# these are small, about 5 cm:
Q1            Larger lunar elliptic diurnal constituent
K2                    Lunisolar semidiurnal constituent
"""
df1['Obs_0'] = df1['Tide Seattle'] - eta0
df3['SNE'] = df3['M2'] + df3['S2'] + df3['N2']
sun_moon_fac = ( df_cons.loc['P1','Amplitude'] /
    (df_cons.loc['P1','Amplitude'] + df_cons.loc['O1','Amplitude']) )
df3['SunDec'] = sun_moon_fac *df3['K1'] + df3['P1']
df3['MoonDec'] = (1-sun_moon_fac)*df3['K1'] + df3['O1']
df3['Trial Sum'] = df3['SNE'] + df3['SunDec'] + df3['MoonDec']

# PLOTTING

if False:
    # initializing
    plt.close('all')
    fig = plt.figure(figsize=(14,8))
    ax = fig.add_subplot(111)
    lw = .5
    df1.plot(ax=ax, lw=lw)
    df3.plot(y='Tide Prediction', ax=ax, lw=lw)
    plt.show()
    
if True:
    # initializing
    plt.close('all')
    fig = plt.figure(figsize=(14,8))
    ax = fig.add_subplot(111)
    lw = .5
    df1.plot(y = 'Obs_0', ax=ax, lw=lw)
    df3.plot(y='SNE', ax=ax, lw=lw)
    df3.plot(y='SunDec', ax=ax, lw=lw)
    df3.plot(y='MoonDec', ax=ax, lw=lw)
    #df3.plot(y='Trial Sum', ax=ax, lw=lw)
    plt.show()
