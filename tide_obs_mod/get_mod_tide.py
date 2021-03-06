"""
Code to automate getting year-long tide height records from
an extracted model NetCDF file of surface height at
a series of NOAA and DFO sites around the Salish Sea and NE Pacific
coast.

Also computes harmonics using utide.

Have to run zeta_extractor.py first.

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import netCDF4 as nc

from importlib import reload
import obsfun as ofn
reload(ofn)

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# standard arguments
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='cas4')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='v2')
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default='lo6biom')
parser.add_argument('-0', '--date_string0', nargs='?', type=str, default='2017.07.20')
parser.add_argument('-1', '--date_string1', nargs='?', type=str, default='2017.07.22')
args = parser.parse_args()

# save some arguments
Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name
ds0 = args.date_string0
ds1 = args.date_string1
Ldir['date_string0'] = ds0
Ldir['date_string1'] = ds1

import zfun

pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun


dir0 = Ldir['parent'] + 'ptools_output/tide/'

noaa_sn_dict, dfo_sn_dict, sn_dict = ofn.get_sn_dicts()
testing = False
if testing == True:
    name_list = ['Seattle']
    a = dict()
    for name in name_list:
        a[name] = sn_dict[name]
    sn_dict = a

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
    hfn = obs_dir + 'h_' + str(sn) + '_' + str(year) + '.p'
    Tobs[name] = pd.read_pickle(fn)
    Mobs[name] = Lfun.csv_to_dict(mfn)
    Hobs[name] = pickle.load(open(hfn, 'rb'))

def get_ij_good(lon, lat, xvec, yvec, i0, j0):
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
    dddm = ddd[mm] # vector of good distances
    # indices of best point
    igood = ii[mm][dddm==dddm.min()][0]
    jgood = jj[mm][dddm==dddm.min()][0]
    #
    return igood, jgood

# set the locations for model input and output
mod_string = Ldir['gtagex'] + '_' + ds0 + '_' + ds1
mod_dir_in = Ldir['LOo'] + 'layer/' + mod_string + '/'

mod_dir_out0 = dir0 + 'mod_data/'
Lfun.make_dir(mod_dir_out0)
mod_dir_out = mod_dir_out0 + Ldir['gtagex'] + '/'
Lfun.make_dir(mod_dir_out)
        
# load model data
fn = mod_dir_in + 'zeta_hourly.nc'
ds = nc.Dataset(fn)
lon = ds['lon_rho'][:]
lat = ds['lat_rho'][:]
mask = ds['mask_rho'][:]
xvec = lon[0,:].flatten()
yvec = lat[:,0].flatten()

Tmod = dict()
Mmod = dict()
Hmod = dict()

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
        i0, j0 = get_ij_good(lon, lat, xvec, yvec, i0, j0)
    # put data into a DataFrame
    zm = ds['zeta'][:,j0,i0].squeeze()
    pair = ds['Pair'][:,j0,i0].squeeze()# units are millibars = 1e-3 * 1e5 Pa = 100 Pa
    tm = ds['ocean_time'][:]
    dtm_list = []
    for t in tm:
        dtm_list.append(Lfun.modtime_to_datetime(t))
    dti = pd.to_datetime(dtm_list)
    dti = dti.tz_localize('UTC')
    df = pd.DataFrame(data={'eta':zm}, index = dti)
    df['Pair (mb)'] = pair
    df.index.name = 'Date'
    #
    Tmod[name] = df
    m_dict = dict()
    m_dict['lon'] = xvec[i0]
    m_dict['lat'] = yvec[j0]
    Mmod[name] = m_dict
    Hmod[name]= ofn.get_harmonics(df, float(Mmod[name]['lat']))

# save results to disk, exactly like what we did for the observations
for name in sn_dict.keys():
    sn = sn_dict[name]
    fn = mod_dir_out + 'tide_' + str(sn) + '_' + str(year) + '.p'
    mfn = mod_dir_out + 'm_' + str(sn) + '_' + str(year) + '.csv'
    hfn = mod_dir_out + 'h_' + str(sn) + '_' + str(year) + '.p'
    Tmod[name].to_pickle(fn)
    Lfun.dict_to_csv(Mmod[name], mfn)
    pickle.dump(Hmod[name], open(hfn, 'wb'))
    
if False:
    plt.close('all')
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
        xm = Mmod[name]['lon']
        ym = Mmod[name]['lat']
        ax.plot(xm, ym, '*b')
    #
    plt.show()
    
