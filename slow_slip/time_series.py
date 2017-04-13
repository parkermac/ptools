"""
Code to plot time series from the bottom pressure fields
that I created from tidally-averaged LiveOcean runs.

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart(gridname='cascadia1', tag='base')
import zrfun
import zfun

# this allows me to plot things on fjord without using an X term.
# import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
from datetime import datetime, timedelta
import pickle
import netCDF4 as nc

import matplotlib.dates as mdates

pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

# specify input
Ldir['gtagex'] = Ldir['gtag'] + '_lobio1'

yr = 2015
out_tag = 'full_' + str(yr)

# setup output location
in_dir0 = Ldir['parent'] + '/ptools_output/slow_slip/'
in_dir = in_dir0 + Ldir['gtagex'] + '_' + out_tag + '/'

bp_arr = pickle.load(open(in_dir + 'bottom_pressure_array.p', 'rb'))
G = pickle.load(open(in_dir + 'grid_info_dict.p', 'rb'))
bp_mean = np.mean(bp_arr, axis=0)
bp_anom = bp_arr - bp_mean

if yr == 2013:
    d0 = 2
else:
    d0 = 1
dt0 = datetime(yr,1,d0,12,0,0)
dt_list = []
for d in range(bp_arr.shape[0]):
    dt = dt0 + timedelta(days = d)
    dt_list.append(dt)
dt_vec = mdates.date2num(dt_list)


if True:
    # find the indices of a given depth contour at many latitudes
    h = G['h']
    x = G['lon_rho']
    y = G['lat_rho']
    ii_list = []
    jj_list = []
    h_target = 500
    for jj in range(0,G['M'],20):
        jj_list.append(jj)
        ii = zfun.find_nearest_ind(h[jj,:], h_target)
        ii_list.append(ii)
    pp = np.nan * np.ones((len(dt_vec), len(ii_list)))
    hh_list = []
    for aa in range(len(jj_list)):
        pp[:,aa] = bp_anom[:, jj_list[aa], ii_list[aa]]
        hh_list.append(h[jj_list[aa], ii_list[aa]])
    hh = np.array(hh_list)
else:
    # or work on a singla latitude
    jj = zfun.find_nearest_ind(G['lat_rho'][:,0],47)
    nstep = 2
    pp = bp_anom[:,jj,::nstep]
    hh = G['h'][jj,::nstep]

# plotting

plt.close()

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)

h_list = []
flag = True
P0 = pp[:,10]
h0 = hh[10]
for ii in range(pp.shape[1]):
    P = pp[:, ii]
    #if P.mask.all() == False and (hh[ii] < 1500) and (hh[ii] > 250):
    #if (hh[ii] < 1500) and (hh[ii] > 250):
    # if flag == True:
    #     P0 = P.copy()
    #     h0 = hh[ii]
    #     flag = False
    ax.plot(dt_list, P-P0, '-')
    h_list.append(str(int(hh[ii])))

ax.set_xlabel('Date')
ax.set_ylabel('Pressure [Pa] (100 Pa = 1 cm)')
ax.set_title('Botom Pressure Anomalies Relative to ' + str(int(h0)))
plt.legend(h_list)

plt.show()

