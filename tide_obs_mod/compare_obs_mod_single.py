"""
Compare detailed tide records at a selected station.

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

name = 'Victoria Harbour'
year  = 2017
dt0 = datetime(year,1,1)
dt1 = datetime(year,1,23)

# observed tide
obs_dir = dir0 + 'obs_data/'
sn = sn_dict[name]
fn = obs_dir + 'tide_' + str(sn) + '_' + str(year) + '.p'
Tobs = pd.read_pickle(fn)
mfn = obs_dir + 'm_' + str(sn) + '_' + str(year) + '.csv'
Mobs = Lfun.csv_to_dict(mfn)

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
def get_model_records(sn_dict, Mobs):
    
    # mod_dir = dir0 + 'mod_data/'
    # fn = mod_dir + gtagex + '/eta_' + str(year) + '.nc'
    
    # hack
    fn = '/Users/pm7/Documents/LiveOcean_output/layer/cas6_v1_lo8_2017.01.01_2017.01.23/zeta_hourly.nc'
    #fn = '/Users/pm7/Documents/LiveOcean_output/layer/cas4_v2_lo6biom_2017.01.01_2017.12.31/zeta_hourly.nc'
    #fn = '/Users/pm7/Documents/LiveOcean_output/layer/cas5_v3_lo8_2017.01.01_2017.12.31/zeta_hourly.nc'
    
    ds = nc.Dataset(fn)
    lon = ds['lon_rho'][:]
    lat = ds['lat_rho'][:]
    mask = ds['mask_rho'][:]
    xvec = lon[0,:].flatten()
    yvec = lat[:,0].flatten()


    slon = Mobs['lon']
    slat = Mobs['lat']
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
    Tmod = df
    m_dict = dict()
    m_dict['lon'] = xvec[i0]
    m_dict['lat'] = yvec[j0]
    Mmod = m_dict
        
    ds.close()
    return Tmod, Mmod
    
Tmod, Mmod = get_model_records(sn_dict, Mobs)

# PLOTTING


Tobs = Tobs[dt0:dt1]
Tmod = Tmod[dt0:dt1]

Tobs = Tobs.rename(columns={'eta':'obs'})
Tmod = Tmod.rename(columns={'eta':'mod'})

tcomb = Tobs.copy()
tcomb['obs'] -= tcomb.loc[:, 'obs'].mean()

tcomb['mod'] = Tmod['mod']
    
tcomb['mod'] -= tcomb.loc[:, 'mod'].mean()
        
tcomb.plot(title=name)
    
plt.show()
    
