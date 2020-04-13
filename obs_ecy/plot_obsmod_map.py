"""
Plots data from the Department of Ecology, combining data from CTD casts
and bottles at the same station.  Also including model output.

This makes a summary map.

"""
import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
import zfun

sys.path.append(os.path.abspath('../../LiveOcean/plotting'))
import pfun

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seawater as sw

Ldir = Lfun.Lstart()
Ldir['gtagex'] = 'cas6_v3_lo8b'
year = 2017
year_str = str(year)

# input
dir11 = Ldir['parent'] + 'ptools_output/ecology/'
out_fn = 'ObsMod_' + Ldir['gtagex'] + '_'+ str(year) + '.p'
Bc = pd.read_pickle(dir11 + out_fn)

# Bc is a DataFrame with processed bottles and casts, both observed and modeled
# at specific depths:
#      Station       Date  Znom    Z  ... Mod DO (mg L-1) Mod DIN (uM) Density (kg m-3) Mod Density (kg m-3)
# 0     ADM001 2017-02-23     0  NaN  ...         10.0564      16.9518          1022.36              1021.47
# 1     ADM001 2017-02-23   -10  NaN  ...          9.2235      20.4011          1022.59              1022.89
# 2     ADM001 2017-02-23   -30  NaN  ...         8.91924      21.2567          1022.67              1023.31
# 3     ADM001 2017-02-23   -80  NaN  ...         8.79093      20.8736          1023.34              1023.72
# The full ist of columns is:
# Index(['Station', 'Date', 'Znom', 'Z', 'Salinity', 'Temp. (deg C)',
#        'Chl (mg m-3)', 'DO (mg L-1)', 'DIN (uM)', 'Mod Salinity',
#        'Mod Temp. (deg C)', 'Mod Chl (mg m-3)', 'Mod DO (mg L-1)',
#        'Mod DIN (uM)', 'Density (kg m-3)', 'Mod Density (kg m-3)'],
#       dtype='object')

Bc.index = Bc.Date

dir0 = '../../ptools_data/ecology/'
# load processed station info and data
sta_df = pd.read_pickle(dir0 + 'sta_df.p')
# add Canadian data
dir1 = '../../ptools_data/canada/'
# load processed station info and data
sta_df_ca = pd.read_pickle(dir1 + 'sta_df.p')
# merge
sta_df = pd.concat((sta_df, sta_df_ca), sort=False)

# PLOTTING
plt.close('all')
fig = plt.figure(figsize=(20,12))
fs = 18
co = 'purple'
cm = 'magenta'
abc = 'abcd'

# axis limits
x0 = -124.5; x1 = -122; y0 = 46; y1 = 49.5 # Salish Sea
aa = [x0, x1, y0, y1]

def msz(vv, v0, v1, scl=40):
    msz = scl*(vv - v0)/(v1 - v0)
    if msz < 1:
        msz = 1
    return msz

#vname = 'Salinity'; v0 = 24; v1 = 34
#vname = 'Temp. (deg C)'; v0 = 8; v1 = 18
#vname = 'DO (mg L-1)'; v0 = 0; v1 = 12
vname = 'DIN (uM)'; v0 = 0; v1 = 35
counter = 1
for z in [0, -30, -30]:
    
    ax = fig.add_subplot(1, 3, counter)
    pfun.add_coast(ax, color='gray')
    ax.axis(aa)
    pfun.dar(ax)
    ax.tick_params(labelsize=.8*fs)
    
    if counter > 1:
        ax.set_yticklabels([])

    for sn in sta_df.index:
        xs = sta_df.loc[sn, 'Longitude']
        ys = sta_df.loc[sn, 'Latitude']
    
    
        mo0 = 6; mo1 = 9 # month range to consider
        mask = ((Bc.Station==sn) & (Bc.Znom==z)
            & (Bc.index.month >= mo0) & (Bc.index.month <= mo1))
    
        vo = Bc.loc[mask,vname].mean()
        vm = Bc.loc[mask,'Mod ' + vname].mean()
        # print('%s vo=%0.2f vm=%0.2f' % (sn, vo, vm))
    
        ax.plot(xs, ys, 'o', markerfacecolor='None', markeredgecolor=co,
            markersize=msz(vo, v0, v1), markeredgewidth=2)
        ax.plot(xs, ys, 'o', markerfacecolor='None', markeredgecolor=cm,
            markersize=msz(vm, v0, v1), markeredgewidth=2)
        ax.text(xs, ys+.05, sn, ha='center', va='center', size=.5*fs)
        
    ax.text(.95, .95, '(%s) Z = %d (m)' % (abc[counter-1], z),
        ha='right', va='center', size=fs, transform=ax.transAxes, weight='bold')
    
    
    if counter == 1:
        # legend
        ax.text(.95, .03, '%s\nMonths %d to %d' % (vname, mo0, mo1),
            ha='right', va='bottom', size=1.3*fs, transform=ax.transAxes)
        ax.plot(-122.5, 46.8, 'o', markerfacecolor='None', markeredgecolor='b',
            markersize=msz(v0, v0, v1), markeredgewidth=2)
        ax.plot(-122.5, 46.8, 'o', markerfacecolor='None', markeredgecolor='b',
            markersize=msz(v1, v0, v1), markeredgewidth=2)
        ax.text(-122.5, 46.6, 'Range = %d to %d' % (v0, v1),
            ha='center', va='center', size=.8*fs, style='italic', color='b')
        ax.text(.05, .06, 'Observed', color=co,
            size=fs, transform=ax.transAxes)
        ax.text(.05, .03, 'Model', color=cm,
            size=fs, transform=ax.transAxes)
            
    counter += 1
    
    
#fig.tight_layout()
#plt.savefig(out_fn)
plt.show()



