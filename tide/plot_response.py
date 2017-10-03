"""
Code to plot the observed and modeled harmonic constituent
amplitude and phase.

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import zfun
import zrfun

pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np
import pickle
import netCDF4 as nc

from importlib import reload
import obsfun as ofn
reload(ofn)

home = os.environ.get('HOME')
dir00 = home + '/Documents/'
dir0 = dir00 + 'ptools_output/tide/'

noaa_sn_dict, dfo_sn_dict, sn_dict = ofn.get_sn_dicts()
testing = False
if testing == True:
    name_list = ['Neah Bay', 'Campbell River', 'Victoria Harbour', 'Tacoma']
    a = dict()
    for name in name_list:
        a[name] = sn_dict[name]
    sn_dict = a

# load observational data
year  = 2013

obs_dir = dir0 + 'obs_data/'
Tobs = dict()
Mobs = dict()
Hobs = dict()
for name in sn_dict.keys():
    # observed tide
    sn = sn_dict[name]
    fn = obs_dir + 'tide_' + str(sn) + '_' + str(year) + '.p'
    mfn = obs_dir + 'm_' + str(sn) + '_' + str(year) + '.csv'
    hfn = obs_dir + 'h_' + str(sn) + '_' + str(year) + '.p'
    #Tobs[name] = pd.read_pickle(fn)
    Mobs[name] = Lfun.csv_to_dict(mfn)
    Hobs[name] = pickle.load(open(hfn, 'rb'))
    
mod_dir = dir0 + 'mod_data/'
Tmod = dict()
Mmod = dict()
Hmod = dict()
for name in sn_dict.keys():
    # observed tide
    sn = sn_dict[name]
    fn = mod_dir + 'tide_' + str(sn) + '_' + str(year) + '.p'
    mfn = mod_dir + 'm_' + str(sn) + '_' + str(year) + '.csv'
    hfn = mod_dir + 'h_' + str(sn) + '_' + str(year) + '.p'
    #Tmod[name] = pd.read_pickle(fn)
    Mmod[name] = Lfun.csv_to_dict(mfn)
    Hmod[name] = pickle.load(open(hfn, 'rb'))
    
# # the main diurnals:
# # - can make into SunDec and MoonDec
# O1                            Lunar diurnal constituent
# P1                            Solar diurnal constituent
# K1                            Lunar diurnal constituent (really Luni-Solar)
#
# # the main semi-diurnals:
# # - can make into Spring-Neap, with elliptical modulation
# M2              Principal lunar semidiurnal constituent
# S2              Principal solar semidiurnal constituent
# N2        Larger lunar elliptic semidiurnal constituent
#
# # these are small, about 5 cm:
# Q1            Larger lunar elliptic diurnal constituent
# K2                    Lunisolar semidiurnal constituent

def get_AG(name, hn, Hobs, Hmod):
    ho = Hobs[name]
    hm = Hmod[name]
    Ao = ho.A[ho.name==hn]
    Am = hm.A[hm.name==hn]
    Go = ho.g[ho.name==hn]
    Gm = hm.g[hm.name==hn]
    Fo = 24*ho.aux.frq[ho.name==hn] # cycles per day
    Fm = 24*hm.aux.frq[hm.name==hn]
    #
    return Ao, Am, Go, Gm, Fo, Fm
    
# plotting

plt.close('all')

if True:
    
    #name0 = 'Westport'
    #name0 = 'La Push'
    name0 = 'Neah Bay'
    #name0 = 'Bamfield'
    
    #name1 = 'Campbell River'
    #name1 = 'Point Atkinson'
    #name1 = 'Vancouver'
    name1 = 'Tacoma'
    #name1 = 'Seattle'
    
    #
    fig = plt.figure(figsize=(8, 8))
    ax1 = fig.add_subplot(211)
    ax1.set_xlim(0, 3)
    ax1.set_ylim(0, 3)
    ax1.grid()
    ax1.set_ylabel('Amplification Factor')
    ax1.set_title(name0 + ' to ' + name1)
    
    ax2 = fig.add_subplot(212)
    ax2.set_xlim(0, 3)
    ax2.set_ylim(0, 1)
    ax2.grid()
    ax2.set_xlabel('Frequency (cycles/day)')
    ax2.set_ylabel('Phase shift / 180 (deg/deg)')
    #
    hn_list = ['M2','S2','N2','O1','P1','K1']
    for hn in hn_list:
        Ao0, Am0, Go0, Gm0, Fo0, Fm0 = get_AG(name0, hn, Hobs, Hmod)
        Ao1, Am1, Go1, Gm1, Fo1, Fm1 = get_AG(name1, hn, Hobs, Hmod)
        Aro = Ao1/Ao0
        Arm = Am1/Am0
        dGo = Go1 - Go0
        if dGo < 0:
            dGo += 360
        dGm = Gm1 - Gm0
        if dGm < 0:
            dGm += 360
        ax1.plot(Fo0, Aro, '*r', Fm0, Arm, '*b')
        ax1.plot(Fo0, Am0/Ao0, 'ok', alpha=.3)
        ax2.plot(Fo0, dGo/180, '*r', Fm0, dGm/180, '*b')
        #
    
if False:
    #
    fig2 = plt.figure(figsize=(8, 8))
    ax = fig2.add_subplot(111)
    pfun.add_coast(ax)
    ax.set_xlim(-130, -122)
    ax.set_ylim(42, 52)
    pfun.dar(ax)
    #
    for name in sn_dict.keys():
        xo = Mobs[name]['lon']
        yo = Mobs[name]['lat']
        ax.plot(xo, yo, '*r')
        ax.text(xo, yo, name)
        xm = Mmod[name]['lon']
        ym = Mmod[name]['lat']
        ax.plot(xm, ym, '*b')
    
plt.show()