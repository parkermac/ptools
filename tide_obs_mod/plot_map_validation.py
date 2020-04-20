"""
Compare observed and modeled amplitude and phase for two constituents
at all stations.  The goal is to create maps that conveys the tidal
validation results.

"""

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart()

sys.path.append(os.path.abspath('../../LiveOcean/plotting'))
import pfun

import pandas as pd
import pickle
import numpy as np
import matplotlib.pyplot as plt

import obsfun as ofn

# select model run, year, and constituents to plot (2)
gtagex = 'cas6_v3_lo8b'
year  = 2017
const_list = ['M2', 'K1']

# input
dir0 = Ldir['parent'] + 'ptools_output/tide/'
obs_dir = dir0 + 'obs_data/'
mod_dir = dir0 + 'mod_data/' + gtagex + '/'

# output
outdir = dir0 + 'obs_mod_summary_output_cas6/'
Lfun.make_dir(outdir)
outfn = outdir + 'map_validation_' + gtagex + '_' + str(year) + '.png'

# get info
noaa_sn_dict, dfo_sn_dict, sn_dict = ofn.get_sn_dicts()
sn_list = sn_dict.keys()

map_df_dict = dict()
for const in const_list:
    map_df = pd.DataFrame(index = sn_list)
    for name in sn_list:
        # load observational data
        sn = sn_dict[name] # station number
        hfn = obs_dir + 'h_' + str(sn) + '_' + str(year) + '.p'
        Hobs = pickle.load(open(hfn, 'rb'))
        # get station locations
        mfn = obs_dir + 'm_' + str(sn) + '_' + str(year) + '.csv'
        M = Lfun.csv_to_dict(mfn) # has name and lon,lat info for a station
        map_df.loc[name,'lon'] = float(M['lon'])
        map_df.loc[name,'lat'] = float(M['lat'])
        # load model data
        hfn = mod_dir + 'h_' + str(sn) + '_' + str(year) + '.p'
        Hmod = pickle.load(open(hfn, 'rb'))
        # get constituent info
        for hn in ofn.hn_list: # hn is a constituent name like M2
            Ao, Am, Go, Gm, Fo, Fm = ofn.get_AG(hn, Hobs, Hmod)
            if hn == const:
                map_df.loc[name,'Ao'] = Ao
                map_df.loc[name,'Am'] = Am
                map_df.loc[name,'Go'] = Go
                map_df.loc[name,'Gm'] = Gm
                map_df.loc[name,'Dph'] = Gm - Go
    map_df_dict[const] = map_df

# plotting
plt.close('all')
fig = plt.figure(figsize=(8,8))
fs = 18
abc = 'abc'
scl = 50

# axis limits
x0 = -127; x1 = -121.5; y0 = 42.5; y1 = 50.5
dx = x1 - x0; dy = y1 - y0

co = 'c'
cm = 'b'

counter = 1
for const in const_list:
    map_df = map_df_dict[const]
    ax = fig.add_subplot(1,2,counter)
    pfun.add_coast(ax, color='gray')
    ax.axis([x0, x1, y0, y1])
    pfun.dar(ax)
    for name in map_df.index:
        x = map_df.loc[name,'lon']
        y = map_df.loc[name,'lat']
        Ao = map_df.loc[name,'Ao']
        Am = map_df.loc[name,'Am']
        Go = map_df.loc[name,'Go']
        Gm = map_df.loc[name,'Gm']
        Dph = map_df.loc[name,'Dph']
        ax.plot(x,y,'o', markersize=Ao*scl,
            markerfacecolor='None', markeredgecolor=co, markeredgewidth=3)
        ax.plot(x,y,'o', markersize=Am*scl,
            markerfacecolor='None', markeredgecolor=cm, markeredgewidth=1)
        # fractional positions
        xx = (x - x0)/dx
        yy = (y - y0)/dy
        ax.quiver(xx,yy, 0*xx, np.ones_like(yy),
            transform=ax.transAxes, scale=20, scale_units='height',
            headwidth=0,headlength=0, angles=Go, color=co)
        ax.quiver(xx,yy, 0*xx, np.ones_like(yy),
            transform=ax.transAxes, scale=20, scale_units='height',
            headwidth=0,headlength=0, angles=Gm, color=cm)
    ax.text(.96, .93, '(%s) $%s_{%s}$' % (abc[counter-1], const[0], const[1]),
        fontsize=fs*1.5, transform=ax.transAxes, ha='right')
    if counter == 1:
        ax.set_xlabel('Longitude', fontsize=fs)
        ax.set_ylabel('Latitude', fontsize=fs)
        ax.text(.03, .08, 'Observed',
            fontsize=fs, transform=ax.transAxes, color=co, style='italic')
        ax.text(.03, .04, 'Modeled',
            fontsize=fs, transform=ax.transAxes, color=cm, style='italic')
    if counter == 2:
        ax.set_xlabel('Longitude', fontsize=fs)
        ax.set_yticklabels([])
        
        # scale
        xx0 = .85; yy0=.2
        ax.plot(xx0,yy0,'o', markersize=scl,
            markerfacecolor='None', markeredgecolor='k',transform=ax.transAxes)
        ax.text(xx0,yy0+.08,'$1\ m$\nAmplitude', size=.6*fs, ha='center', va='center',transform=ax.transAxes)
        ax.quiver(xx0,yy0, 0, 1,
            transform=ax.transAxes, scale=20, scale_units='height',
            headwidth=0,headlength=0, angles=0, color=co)
        ax.quiver(xx0,yy0, 0, 1,
            transform=ax.transAxes, scale=20, scale_units='height',
            headwidth=0,headlength=0, angles=15, color=cm)
        ax.text(xx0,yy0-.08,'Model Phase', size=.6*fs, ha='center', va='center',transform=ax.transAxes, color=cm)
        ax.text(xx0,yy0-.11,'lags by $15^{\circ}$', size=.6*fs, ha='center', va='center',transform=ax.transAxes)
        
    ax.tick_params(labelsize=.8*fs)
    counter += 1
    
fig.tight_layout()
plt.savefig(outfn)
plt.show()
