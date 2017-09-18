"""
Plot of Decadal changes in predicted tide.
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

indir = os.environ.get('HOME') + '/Documents/ptools_data/tide/'
city = 'Seattle'

def make_tide(year, city):
    # MAKE A SYNTHETIC TIDE

    # and set related time limits
    dt0 = datetime(year,1,1,tzinfo=pytz.timezone('UTC'))
    dt1 = datetime(year,12,31,tzinfo=pytz.timezone('UTC'))
    index = pd.date_range(dt0, dt1, freq='h')

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
    cons_fn = indir + city.lower() + '_constituents.txt'
    cons_df = pd.read_csv(cons_fn, header=0, sep='\t', index_col='Name')
    #
    # Create the synthetic tide (inherits tzinfo from obs_df)
    pred_df = pd.DataFrame(index=index, columns=cons_list)
    pred_df['tsec'] = (pred_df.index - dt0).total_seconds()
    # separate constituents
    for cons in cons_list:
        f = f_ser[cons]
        H = cons_df.loc[cons,'Amplitude']
        om = (np.pi/180) * cons_df.loc[cons, 'Speed']/3600
        vu = (np.pi/180) * vu_ser[cons]
        G = (np.pi/180) * cons_df.loc[cons, 'Phase']
        this_eta = f * H * np.cos(om * pred_df['tsec'] + vu - G)
        pred_df[cons] = this_eta * m2f
    # create full 8-constituent prediction
    pred_df['z'] = 0
    for cons in cons_list:
        pred_df['z'] += pred_df[cons]
    # add the mean se level at Seattle (feet)
    pred_df['z'] += 6.947
    d = np.array(index.dayofyear.tolist())
    h = np.array(index.hour.tolist())
    yd = d + h/24
    pred_df['yd'] = yd
    return pred_df
    
df_2017 = make_tide(2017, city)
df_2025 = make_tide(2025, city)

df_2017.index=df_2017['yd']
df_2025.index=df_2025['yd']

# PLOTTING
plt.close('all')

# RC SETUP (plotting defaults)
def set_rc(fs, lw, mks):
    fs_big = fs
    fs_small = fs-4
    lw_big = lw
    lw_small = lw+1
    plt.rc('xtick', labelsize=fs_small)
    plt.rc('ytick', labelsize=fs_small)
    plt.rc('xtick.major', size=10, pad=5, width=lw_small)
    plt.rc('ytick.major', size=10, pad=5, width=lw_small)
    plt.rc('axes', lw=lw_small)
    plt.rc('lines', lw=lw_big, markersize=mks)
    plt.rc('font', size=fs_big)
    plt.rc('grid', color='g', ls='-', lw=lw_small, alpha=.3)
fs = 16
lw = 3
lwt = 2
mks = 25
set_rc(fs, lw, mks)

fig = plt.figure(figsize=(18,8))
ax = fig.add_subplot(111)
df_2025.plot(y='z',
        legend=False, style='-r', ax=ax, lw=lwt, grid=True)
df_2017.plot(y='z', title=('Predicted Tide Height (feet) ' + city),
        legend=False, style='-b', ax=ax, lw=lwt, grid=True, alpha=.5)
        
ax.text(.23,.11,'2017', color='b',
    transform=ax.transAxes, fontweight='bold', fontsize=26)
ax.text(.27,.03,'2025', color='r',
    transform=ax.transAxes, fontweight='bold', fontsize=26)

ax.set_xlabel('Day of Year')    
plt.show()

# RC CLEANUP
plt.rcdefaults()
