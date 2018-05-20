"""
Compare detailed tide records at a selected point

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

from importlib import reload
import obsfun as ofn
reload(ofn)

home = os.environ.get('HOME')
dir00 = home + '/Documents/'
dir0 = dir00 + 'ptools_output/tide/'

noaa_sn_dict, dfo_sn_dict, sn_dict = ofn.get_sn_dicts()

# Station pairs, North-to-South:
#
# - Strait of Georgia
#name_list = ['La Push', 'Campbell River']
#name_list = ['La Push', 'Point Atkinson']
name_list = ['La Push', 'Vancouver']
# - JdF
#name_list = ['La Push', 'Victoria Harbour']
# - Puget Sound
#name_list = ['La Push', 'Seattle']
#name_list = ['La Push', 'Tacoma']

a = dict()
for name in name_list:
    a[name] = sn_dict[name]
sn_dict = a

# select several model runs
run_list = ['cas4_v0_lo6m_goodtide']

# load observational data
year  = 2017
obs_dir = dir0 + 'obs_data/'
Tobs = dict()
Mobs = dict()
Hobs = dict()
for name in sn_dict.keys():
    # observed tide
    sn = sn_dict[name]
    fn = obs_dir + 'tide_' + str(sn) + '_' + str(year) + '.p'
    mfn = obs_dir + 'm_' + str(sn) + '_' + str(year) + '.csv'
    Tobs[name] = pd.read_pickle(fn)
    Mobs[name] = Lfun.csv_to_dict(mfn)

def get_ij_good(lon, lat, xvec, yvec, i0, j0, mask):
    # find the nearest unmasked point
    #
    # starting point
    lon0 = lon[j0,i0]
    lat0 = lat[j0,i0]
    pad = 5 # how far to look (points)
    # indices of box to search over
    imax = len(xvec)-1
    jmax = len(yvec)-1
    I = np.arange(i0-pad, i0+pad)
    J = np.arange(j0-pad,j0+pad)
    # account for out-of-range points
    if I[0] < 0:
        I = I - I[0]
    if I[-1] > imax:
        I = I - (I[-1] - imax)
    if J[0] < 0:
        J = J - J[0]
    if J[-1] > jmax:
        J = J - (J[-1] - jmax)
    ii, jj = np.meshgrid(I, J)
    # sub arrays
    llon = lon[jj,ii]
    llat = lat[jj,ii]
    xxx, yyy = zfun.ll2xy(llon, llat, lon0, lat0)
    ddd = np.sqrt(xxx**2 + yyy**2) # distance from original point
    mmask = mask[jj,ii] 
    mm = mmask==1 # Boolean array of good points
    #print(mm)
    dddm = ddd[mm] # vector of good distances
    # indices of best point
    igood = ii[mm][dddm==dddm.min()][0]
    jgood = jj[mm][dddm==dddm.min()][0]
    #
    return igood, jgood

# function to load model data
def get_model_records(gtagex, dir0, year, sn_dict, Mobs):
    
    mod_dir = dir0 + 'mod_data/'
    fn = mod_dir + gtagex + '/eta_' + str(year) + '.nc'
    ds = nc.Dataset(fn)
    lon = ds['lon_rho'][:]
    lat = ds['lat_rho'][:]
    mask = ds['mask_rho'][:]
    xvec = lon[0,:].flatten()
    yvec = lat[:,0].flatten()

    Tmod = dict()
    Mmod = dict()

    for name in sn_dict.keys():
        slon = Mobs[name]['lon']
        slat = Mobs[name]['lat']
        i0, i1, frx = zfun.get_interpolant(np.array([float(slon)]), xvec)
        j0, j1, fry = zfun.get_interpolant(np.array([float(slat)]), yvec)
        i0 = int(i0)
        j0 = int(j0)
        # find indices of nearest good point
        if mask[j0,i0] == 1:
            print(name + ': point ok')
        elif mask[j0,i0] == 0:
            print(name + ':point masked')
            i0, j0 = get_ij_good(lon, lat, xvec, yvec, i0, j0, mask)
        # put data into a DataFrame
        zm = ds['zeta'][:,j0,i0].squeeze()
        tm = ds['ocean_time'][:]
        dtm_list = []
        for t in tm:
            dtm_list.append(Lfun.modtime_to_datetime(t))
        dti = pd.to_datetime(dtm_list)
        dti = dti.tz_localize('UTC')
        df = pd.DataFrame(data={'eta':zm}, index = dti)
        df.index.name = 'Date'
        #
        Tmod[name] = df
        m_dict = dict()
        m_dict['lon'] = xvec[i0]
        m_dict['lat'] = yvec[j0]
        Mmod[name] = m_dict
    ds.close()
    return Tmod, Mmod
    
# get the model data
Tmod_dict = dict()
Mmod_dict = dict()
for gtagex in run_list:
    Tmod_dict[gtagex], Mmod_dict[gtagex] = get_model_records(
        gtagex, dir0, year, sn_dict, Mobs)

# PLOTTING
plt.close('all')

fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True,
                     figsize=(12,8), squeeze=False)
# NOTE: for some reason, when I use sharex=True and
# have a Canadian station FIRST in the name_list the
# first plot is blank - even though the data is good.
# Maybe a problem with the time index?
ax1 = axes[0,0]
ax2 = axes[1,0]

dt0 = datetime(year,1,1)
dt1 = datetime(year,1,31)

count = 0
for name in name_list:#sn_dict.keys():
    tobs = Tobs[name]
    tcomb = tobs.copy()
    tcomb = tcomb.rename(columns={'eta':'eta_obs'})

    eraw = tcomb.loc[dt0:dt1, 'eta_obs'].values
    efilt = zfun.filt_godin(eraw)
    enew = eraw - efilt
    tcomb.loc[dt0:dt1, 'eta_obs'] = enew
    
    #tcomb['eta_obs'] -= tcomb.loc[dt0:dt1, 'eta_obs'].mean()
    
    for gtagex in run_list:
        Tmod = Tmod_dict[gtagex]
        tmod = Tmod[name]
        tcomb[gtagex] = tmod['eta']
        #tcomb[gtagex] -= tcomb.loc[dt0:dt1, gtagex].mean()
        eraw = tcomb.loc[dt0:dt1, gtagex].values
        efilt = zfun.filt_godin(eraw)
        enew = eraw - efilt
        tcomb.loc[dt0:dt1, gtagex] = enew
    if count == 0:
        ax = ax1
    else:
        ax = ax2
        
    tcomb.plot(ax = ax, title=name)
    ax.set_xlim(dt0, dt1)
    ax.set_ylim(-3,2)
    ax.grid()
    count += 1
    
plt.show()

    
if True:
    #
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
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
        Mmod = Mmod_dict[run_list[0]]
        xm = Mmod[name]['lon']
        ym = Mmod[name]['lat']
        ax.plot(xm, ym, '*b')
    #
    plt.show()
    
