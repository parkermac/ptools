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
import matplotlib.dates as mdates
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

#==============================================
testing = False
for_web = True # plots styled for the validation website

if testing==True:
    sn_list = ['Seattle', 'Tacoma']
    save_fig=False

else:
    sn_list = sn_dict.keys()
    save_fig=True
    if for_web==True:
        tag = 'web'
    else:
        tag = 'val'

# where to put plots
if save_fig==True:
    outdir = dir0 + tag + '_series_' + gtagex + '_'+ str(year) + '/'
    Lfun.make_dir(outdir, clean=True)
#==============================================
    
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

# prepare to gather average statistics
err_std_list = []
err_pct_list = []
obs_std_list = []

sn_lonlat_dict = dict()

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
    
    if for_web==True:
        axtup = (2,1)
        fig = plt.figure(figsize=(6,4))
    else:
        axtup = (2,4)
        fig = plt.figure(figsize=(18,8))

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
    tcomb['SSH mod low-passed with Atm. Pressure Effect'] = zfun.filt_godin(tcomb.loc[dt0:dt1, 'SSH mod adj'].values)

    # Hourly SSH for One Month
    ax = plt.subplot2grid(axtup, (0,0), colspan=2)
    if for_web==True:
        tstr = name + ' Tides for one month'
    else:
        tstr = name + ' Tides for one month of ' + str(year) + ' (means removed)'
        
    if for_web == True:
        tcomb.plot(y=['SSH obs', 'SSH mod'], ax = ax, color=['r','b'], x_compat=True)
    else:
        tcomb.plot(y=['SSH obs', 'SSH mod', 'SSH error'], ax = ax, color=['r','b','gray'], x_compat=True)
    ax.legend(loc='lower right', ncol=2)
    ax.set_xlim(dt00, dt11)
    ax.set_ylim(-3.5,3.5)
    ax.set_ylabel('SSH (m)')
    if for_web==True:
        ax.text(.5,.9, tstr, weight='bold', color='k', transform=ax.transAxes, horizontalalignment='center')
    else:
        ax.text(.5,.9, tstr, weight='bold', color='k', transform=ax.transAxes, horizontalalignment='center')
        ax.text(.5, .8, ('Obs Mean = %0.2f (m) Mod Mean = %0.2f (m)' % (obs_mean, mod_mean)),
            transform=ax.transAxes, horizontalalignment='center')
    ax.text(.01, .1, ('SSH Std. Dev. = %0.2f (m)' % (obs_std)),
            transform=ax.transAxes, weight='bold')
    ax.text(.01, .03, ('Error Std. Dev. = %0.2f (m) (%d%%)' % (err, int(100*err/obs_std))),
            transform=ax.transAxes, weight='bold')
    ax.grid(False)
    if for_web==True:
        ax.set_xticklabels([''])
    ax.set_xlabel('')
    
    # save statistics
    err_std_list.append(err)
    err_pct_list.append(100*err/obs_std)
    obs_std_list.append(obs_std)
    
    
    # Full year of low-passed SSH
    ax = plt.subplot2grid(axtup, (1,0), colspan=2)
    tstr = 'Tidally-averaged SSH for full year ' + str(year)
    if for_web==True:
        tcomb.plot(y=['SSH obs low-passed', 'SSH mod low-passed'],
            ax = ax,color=['r','b'], x_compat=True)
    else:
        tcomb.plot(y=['SSH obs low-passed', 'SSH mod low-passed',
            'SSH mod low-passed with Atm. Pressure Effect'],
            ax = ax,color=['r','b','c'], x_compat=True)
    ax.legend(loc='lower right', ncol=2)
    ax.text(.5,.9, tstr, weight='bold', color='k', transform=ax.transAxes, horizontalalignment='center')
    ax.set_xlim(dt0, dt1)
    ax.set_ylim(-.8, .8)
    ax.set_ylabel('SSH (m)')
    
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax.xaxis.set_tick_params(labelrotation=30)
    
    ax.set_xlabel('Date ' + str(year))
    
    
    
    ax.grid()
    
    if for_web == False:
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
            
    # grab location info
    Name = name.replace(' ','_')
    sn_lonlat_dict[Name] = (float(Mobs['lon']), float(Mobs['lat']))

    if for_web==False:
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
    
    if save_fig==False:
        plt.show()
    elif save_fig==True:
        plt.savefig(outdir + Name + '.png')
        plt.close()

# print statistics
print('')
print('Average Error Statistics')
print('RMS Error = %0.2f (m)' % (np.array(err_std_list).mean()))
print('RMS Error = %0.2f (%% of RMS Tide)' % (np.array(err_pct_list).mean()))
print('RMS Tide = %0.2f (m)' % (np.array(obs_std_list).mean()))
print('')

# print list of stations to use with webpage Google Map API
print('sta_list = [')
for station in sn_lonlat_dict.keys():
    lon = sn_lonlat_dict[station][0]
    lat = sn_lonlat_dict[station][1]
    print('    {sta:"%s", lat:%0.4f, lng:%0.4f},' % (station.title(), lat, lon))
print('];')

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
    
