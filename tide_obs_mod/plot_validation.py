"""
Compare model and observations at a single station.

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import zfun

Ldir = Lfun.Lstart()

pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
import netCDF4 as nc
import pickle

from importlib import reload
import obsfun as ofn
reload(ofn)

dir0 = Ldir['parent'] + 'ptools_output/tide/'

# select model run
gtagex = 'cas4_v2_lo6biom'
year  = 2017

noaa_sn_dict, dfo_sn_dict, sn_dict = ofn.get_sn_dicts()

testing = False
if testing==False:
    sn_list = sn_dict.keys()
    save_plot=True
elif testing==True:
    sn_list = ['Seattle']
    save_plot=False

if save_plot==True:
    # place for output
    outdir = dir0 + 'validation_' + gtagex + '_' + str(year) + '/'
    Lfun.make_dir(outdir, clean=True)
    
def get_AG(hn, Hobs, Hmod):
    ho = Hobs
    hm = Hmod
    Ao = ho.A[ho.name==hn]
    Am = hm.A[hm.name==hn]
    Go = ho.g[ho.name==hn]
    Gm = hm.g[hm.name==hn]
    Fo = 24*ho.aux.frq[ho.name==hn] # cycles per day
    Fm = 24*hm.aux.frq[hm.name==hn]
    #
    return Ao, Am, Go, Gm, Fo, Fm

plt.close('all')    
for name in sn_list:

    # load observational data
    obs_dir = dir0 + 'obs_data/'
    sn = sn_dict[name]
    fn = obs_dir + 'tide_' + str(sn) + '_' + str(year) + '.p'
    mfn = obs_dir + 'm_' + str(sn) + '_' + str(year) + '.csv'
    hfn = obs_dir + 'h_' + str(sn) + '_' + str(year) + '.p'
    Tobs = pd.read_pickle(fn)
    Mobs = Lfun.csv_to_dict(mfn)
    Hobs = pickle.load(open(hfn, 'rb'))
        
    # load model data
    mod_dir = dir0 + 'mod_data/' + gtagex + '/'
    sn = sn_dict[name]
    fn = mod_dir + 'tide_' + str(sn) + '_' + str(year) + '.p'
    mfn = mod_dir + 'm_' + str(sn) + '_' + str(year) + '.csv'
    hfn = mod_dir + 'h_' + str(sn) + '_' + str(year) + '.p'
    Tmod = pd.read_pickle(fn)
    Mmod = Lfun.csv_to_dict(mfn)
    Hmod = pickle.load(open(hfn, 'rb'))
    
    # PLOTTING

    fig = plt.figure(figsize=(18,8))
    axtup = (2,4)

    # full year
    dt0 = datetime(year,1,1)
    dt1 = datetime(year,12,31)
    
    # focus times
    dt00 = datetime(year,3,1)
    dt11 = datetime(year,3,31)

    tobs = Tobs
    tcomb = Tobs.copy()
    tcomb = tcomb.rename(columns={'eta':'SSH obs'})
    tcomb['SSH mod'] = Tmod['eta']
    pair = Tmod['Pair (mb)'] - Tmod['Pair (mb)'].mean()
    tcomb['SSH mod adj'] = Tmod['eta'] - pair*100/(1025*9.8)
    tcomb = tcomb.loc[dt0:dt1, :]
    
    # remove means
    obs_mean = tcomb.loc[:, 'SSH obs'].mean()
    mod_mean = tcomb.loc[:, 'SSH mod'].mean()
    moda_mean = tcomb.loc[:, 'SSH mod adj'].mean()
    
    tcomb['SSH obs'] -= obs_mean
    tcomb['SSH mod'] -= mod_mean
    tcomb['SSH mod adj'] -= moda_mean
    tcomb['SSH error'] = tcomb['SSH mod'] - tcomb['SSH obs']
    
    # note that rms is the same as std here because these have ~zero mean
    err = tcomb['SSH error'].std()
    obs_std = tcomb['SSH obs'].std()
    
    tcomb['SSH obs low-passed'] = zfun.filt_godin(tcomb.loc[dt0:dt1, 'SSH obs'].values)
    tcomb['SSH mod low-passed'] = zfun.filt_godin(tcomb.loc[dt0:dt1, 'SSH mod'].values)
    tcomb['SSH mod adj low-passed'] = zfun.filt_godin(tcomb.loc[dt0:dt1, 'SSH mod adj'].values)

    # Hourly SSH for One Month
    ax = plt.subplot2grid(axtup, (0,0), colspan=2)
    tstr = name + ' Tides for one month of ' + str(year) + ' (means removed)'
    tcomb.plot(y=['SSH obs', 'SSH mod', 'SSH error'], ax = ax, color=['r','b','gray'], x_compat=True)
    ax.legend(loc='lower right', ncol=2)
    ax.set_xlim(dt00, dt11)
    ax.set_ylim(-3,3)
    ax.set_ylabel('SSH (m)')
    ax.text(.5,.9, tstr, weight='bold', color='k', transform=ax.transAxes, horizontalalignment='center')
    ax.text(.5, .8, ('Obs Mean = %0.2f (m) Mod Mean = %0.2f (m)' % (obs_mean, mod_mean)),
        transform=ax.transAxes, horizontalalignment='center')
    ax.text(.05, .1, ('SSH Std. Dev. = %0.2f (m)' % (obs_std)),
            transform=ax.transAxes, weight='bold')
    ax.text(.05, .03, ('Error Std. Dev. = %0.2f (m) (%d%%)' % (err, int(100*err/obs_std))),
            transform=ax.transAxes, weight='bold')
    ax.grid(False)
    #ax.set_xticklabels([''])
    ax.set_xlabel('')
    
    # Full year of low-passed SSH
    ax = plt.subplot2grid(axtup, (1,0), colspan=2)
    tstr = 'Tidally-averaged SSH for full year ' + str(year)
    tcomb.plot(y=['SSH obs low-passed', 'SSH mod low-passed','SSH mod adj low-passed'], ax = ax,color=['r','b','c'], x_compat=True)
    ax.legend(loc='lower right', ncol=2)
    ax.text(.5,.9, tstr, weight='bold', color='k', transform=ax.transAxes, horizontalalignment='center')
    ax.set_xlim(dt0, dt1)
    ax.set_ylim(-.8, .8)
    ax.set_ylabel('SSH (m)')
    ax.grid()
    
    # Amplitudes of constituents
    flo = .5
    fhi = 2.5
    ax = plt.subplot2grid(axtup, (0,2), colspan=1)
    ax.set_xlim(flo, fhi)
    ax.set_ylim(0, 1.4)
    ax.grid()
    ax.set_xlabel('Frequency (cycles/day)')
    ax.text(.05,.9, 'Amplitude (m)', weight='bold', color='k',
        transform=ax.transAxes)
    ax.text(.05,.8, 'OBSERVATION', weight='bold', color='r',
        transform=ax.transAxes)
    ax.text(.05,.7, 'MODEL', weight='bold', color='b',
        transform=ax.transAxes)
    #
    hn_list = ['M2','S2','N2','O1','P1','K1']
    for hn in hn_list:
        Ao, Am, Go, Gm, Fo, Fm = get_AG(hn, Hobs, Hmod)
        ax.text(Fo, Ao, hn, color='r', weight='bold',
            horizontalalignment='center', verticalalignment='center')
        ax.text(Fm, Am, hn, color='b', weight='bold',
            horizontalalignment='center', verticalalignment='center')

    # Phases of constituents
    ax = plt.subplot2grid(axtup, (1,2), colspan=1)
    ax.set_xlim(flo, fhi)
    ax.set_ylim(0, 360)
    ax.grid()
    ax.set_xlabel('Frequency (cycles/day)')
    ax.text(.05,.9, 'Phase (degrees)', weight='bold', color='k',
        transform=ax.transAxes)
    #
    hn_list = ['M2','S2','N2','O1','P1','K1']
    for hn in hn_list:
        Ao, Am, Go, Gm, Fo, Fm = get_AG(hn, Hobs, Hmod)
        ax.text(Fo, Go, hn, color='r', weight='bold',
            horizontalalignment='center', verticalalignment='center')
        ax.text(Fm, Gm, hn, color='b', weight='bold',
            horizontalalignment='center', verticalalignment='center')

    # Location Map
    ax = plt.subplot2grid(axtup, (0,3), rowspan=2)
    pfun.add_coast(ax)
    ax.set_xlim(-127, -122)
    ax.set_ylim(43, 51)
    pfun.dar(ax)
    #
    mm = 12
    xo = Mobs['lon']
    yo = Mobs['lat']
    ax.plot(float(xo), float(yo), 'or', markersize=mm)
    xm = Mmod['lon']
    ym = Mmod['lat']
    ax.plot(float(xm), float(ym), '*b', markersize=mm)
    ax.set_title(name)
    
    fig.tight_layout()
    
    if save_plot==False:
        plt.show()
    elif save_plot==True:
        plt.savefig(outdir + name.replace(' ','_') + '.png')
        plt.close()


#
#
# if True:
#     #
#     fig = plt.figure(figsize=(8, 8))
#     ax = fig.add_subplot(111)
#     pfun.add_coast(ax)
#     ax.set_xlim(-130, -122)
#     ax.set_ylim(42, 52)
#     pfun.dar(ax)
#     #
#     for name in sn_dict.keys():
#         xo = Mobs[name]['lon']
#         yo = Mobs[name]['lat']
#         ax.plot(xo, yo, '*r')
#         ax.text(xo, yo, name)
#         Mmod = Mmod_dict[run_list[0]]
#         xm = Mmod[name]['lon']
#         ym = Mmod[name]['lat']
#         ax.plot(xm, ym, '*b')
#     #
#     plt.show()
    
