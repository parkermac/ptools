"""
Compare observed and modeled amplitude and phase for ONE constituents
at all stations.  The goal is to create a single map that conveys the tidal
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

dir0 = Ldir['parent'] + 'ptools_output/tide/'

# select model run
gtagex = 'cas6_v3_lo8b'
year  = 2017

noaa_sn_dict, dfo_sn_dict, sn_dict = ofn.get_sn_dicts()
sn_list = sn_dict.keys()

def get_AG(hn, Hobs, Hmod):
    ho = Hobs
    hm = Hmod
    Ao = ho.A[ho.name==hn]
    Am = hm.A[hm.name==hn]
    Go = ho.g[ho.name==hn]
    Gm = hm.g[hm.name==hn]
    Fo = 24*ho.aux.frq[ho.name==hn] # cycles per day
    Fm = 24*hm.aux.frq[hm.name==hn]
    # Amplitude (A) [m], Phase (G) [deg], Frequency (F) [cpd]
    # o = observed, m = model
    return Ao, Am, Go, Gm, Fo, Fm

sn_coast = ['Charleston', 'South Beach', 'Garibaldi', 'Toke Point',
    'Westport', 'La Push', 'Neah Bay', 'Tofino', 'Bamfield']
sn_salish = ['Port Angeles', 'Friday Harbor', 'Cherry Point', 'Port Townsend',
    'Seattle', 'Tacoma', 'Point Atkinson', 'Vancouver', 'Patricia Bay',
    'Victoria Harbour', 'Campbell River', 'New Westminster']

df_dict = dict() # each DataFrame has one constituent
df = pd.DataFrame(index=sn_list, columns=['ar', 'dph'])
for hn in ofn.hn_list:
    df_dict[hn] = df.copy()

map_df_dict = dict()
const_list = ['M2', 'K1']
#const_list = ['S2', 'O1']
for const in const_list:
    map_df = pd.DataFrame(index = sn_list)#, columns=['lon','lat','Ao','Am','Dph'])
    for name in sn_list:
        # load observational data
        obs_dir = dir0 + 'obs_data/'
        sn = sn_dict[name]
        hfn = obs_dir + 'h_' + str(sn) + '_' + str(year) + '.p'
        Hobs = pickle.load(open(hfn, 'rb'))
        # get station locations
        mfn = obs_dir + 'm_' + str(sn) + '_' + str(year) + '.csv'
        M = Lfun.csv_to_dict(mfn)
        map_df.loc[name,'lon'] = float(M['lon'])
        map_df.loc[name,'lat'] = float(M['lat'])
        # load model data
        mod_dir = dir0 + 'mod_data/' + gtagex + '/'
        hfn = mod_dir + 'h_' + str(sn) + '_' + str(year) + '.p'
        Hmod = pickle.load(open(hfn, 'rb'))
        # get constituent info
        for hn in ofn.hn_list: # hn is a constituent name like M2
            Ao, Am, Go, Gm, Fo, Fm = get_AG(hn, Hobs, Hmod)
            df_dict[hn].loc[name, 'ar'] = Am/Ao
            # fix phase difference when they straddle 360
            if (Gm - Go) > 180:
                Gm = Gm - 360
            elif (Gm - Go) < -180:
                Gm = Gm + 360
            else:
                pass
            df_dict[hn].loc[name, 'dph'] = Gm - Go
            if hn == const:
                map_df.loc[name,'Ao'] = Ao[0]
                map_df.loc[name,'Am'] = Am[0]
                map_df.loc[name,'Go'] = Go[0]
                map_df.loc[name,'Gm'] = Gm[0]
                map_df.loc[name,'Dph'] = Gm[0] - Go[0]
    map_df_dict[const] = map_df

# plotting
plt.close('all')
fig = plt.figure(figsize=(8,8))
fs = 18

x0 = -127; x1 = -121.5
y0 = 42.5; y1 = 50.5
dx = x1 - x0; dy = y1 - y0

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
        scl = 50
        ax.plot(x,y,'o', markersize=Ao*scl,
            markerfacecolor='None', markeredgecolor='k')
        ax.plot(x,y,'o', markersize=Am*scl,
            markerfacecolor='None', markeredgecolor='r')
        # fractional positions
        xx = (x - x0)/dx
        yy = (y - y0)/dy
        ax.quiver(xx,yy, 0*xx, np.ones_like(yy),
            transform=ax.transAxes, scale=20, scale_units='height',
            headwidth=0,headlength=0, angles=Go)
        ax.quiver(xx,yy, 0*xx, np.ones_like(yy),
            transform=ax.transAxes, scale=20, scale_units='height',
            headwidth=0,headlength=0, angles=Gm, color='r')
    ax.text(.05, .9, '$%s_{%s}$' % (const[0], const[1]),
        fontsize=fs*1.5, transform=ax.transAxes)
    if counter == 1:
        ax.set_xlabel('Longitude', fontsize=fs)
        ax.set_ylabel('Latitude', fontsize=fs)
        ax.text(.03, .08, 'Observed',
            fontsize=fs, transform=ax.transAxes, color='k', style='italic')
        ax.text(.03, .04, 'Modeled',
            fontsize=fs, transform=ax.transAxes, color='r', style='italic')
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
            headwidth=0,headlength=0, angles=0)
        ax.quiver(xx0,yy0, 0, 1,
            transform=ax.transAxes, scale=20, scale_units='height',
            headwidth=0,headlength=0, angles=15, color='r')
        ax.text(xx0,yy0-.08,'Model Phase', size=.6*fs, ha='center', va='center',transform=ax.transAxes, color='r')
        ax.text(xx0,yy0-.11,'lags by $15^{\circ}$', size=.6*fs, ha='center', va='center',transform=ax.transAxes)
        
    ax.tick_params(labelsize=.8*fs)
    counter += 1
fig.tight_layout()
plt.show()
