"""
Code to start looking at Slow Slip files.

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()
import zrfun
import matfun
import zfun

import matplotlib.pyplot as plt
import pandas as pd


import numpy as np
from datetime import datetime, timedelta
import pickle
import netCDF4 as nc
import time

pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

#%% setup input locations
R_in_dir0 = Ldir['parent'] + 'roms/output/salish_2006_4_lp/'
# specify ROMS file to work on
dt0 = datetime(2006,1,1)
ndays = 191
# 209 is 2006.07.29, 200 is 2006.07.20
dt = dt0 + timedelta(days=ndays)
f_string = 'f' + dt.strftime('%Y.%m.%d')
print('Working on day ' + f_string)
R_in_dir = R_in_dir0 + f_string + '/'
R_fn = R_in_dir + 'low_passed.nc'
ds = nc.Dataset(R_fn)

in_dir0 = Ldir['parent'] + 'ptools_data/slow_slip/'
in_dir = in_dir0 + 'Fredrickson_Files_2016.11.03/'

sta_fn = in_dir + 'year3stations.mat'
ser_fn = in_dir + 'CIexample.mat'
    
sta = matfun.loadmat(sta_fn)['CIstations']

ser = matfun.loadmat(ser_fn)['exStation']

# 'timeunit': 'second since 1970'
dt_list = []
for t in ser['time']:
    dt_list.append( datetime(1970,1,1) + timedelta(days = t/86400) )

import matplotlib.dates as mdates
mdt = mdates.date2num(dt_list)

days_from_2013 = mdt - mdates.date2num(datetime(2013,1,1))

p = ser['pressure']

pp = ( p - p.mean() ) # "equivalent cm"

# Ser = pd.Series(dict(zip(dt_list,p)))
# Ser.plot()

# load a mooring extraction
m_fn = Ldir['LOo'] + 'moor/cascadia1_base_lobio1_FN01C_low_pass_2013.01.02_2014.12.31.nc'

m_ds = nc.Dataset(m_fn)

ot = m_ds['ocean_time'][:]
m_dt_list = []
for t in ot:
    m_dt_list.append( datetime(1970,1,1) + timedelta(days = t/86400) )
m_mdt = mdates.date2num(m_dt_list)

m_days_from_2013 = m_mdt - mdates.date2num(datetime(2013,1,1))

rho = m_ds['rho'][:] + 1000.
z_w = m_ds['z_w'][:]

g = 9.8
DZ = np.diff(z_w, axis = 0)
bp = (g * DZ * rho).sum(axis = 0) 

m_p = bp - bp.mean()

# convert to equivalent cm
m_p = m_p * 100 / (1025 * g)

#zfun.ncd(M_ds)

plt.close('all')

# PLOT time series

fig = plt.figure(figsize=(10,6))
axs = fig.add_subplot(111)
ppf = zfun.filt_godin(pp)
axs.plot(days_from_2013, zfun.filt_hanning(ppf, n=5*24), '-r', m_days_from_2013, m_p, '-b')
#axs.plot(days_from_2013, pp, '-r', m_days_from_2013, m_p, '-b')

plt.show()

if True:
    # PLOT station map
    aa = [-131, -122, 40, 50]
    fig = plt.figure(figsize=(10,6))
    axm = fig.add_subplot(111)

    axm.plot(sta['lon'], sta['lat'], '*r')

    for ii in range(len(sta['lon'])):
    
        # axm.text(sta['lon'][ii], sta['lat'][ii], sta['name'][ii])
    
        if sta['name'][ii] == 'FN01C':
            axm.text(sta['lon'][ii], sta['lat'][ii], sta['name'][ii])
            print('%s lon = %0.5f lat = %0.5f' % (sta['name'][ii], sta['lon'][ii], sta['lat'][ii]))

    axm.grid()

    axm.axis(aa)

    pfun.add_coast(axm)
    pfun.dar(axm)

    pfun.add_bathy_contours(axm, ds)


    plt.show()

ds.close()