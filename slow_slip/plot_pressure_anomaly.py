"""
Make maps of bottom pressure anomaly.

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart(gridname='cascadia1', tag='base')
import zrfun
import matfun
import zfun
import zrfun

import matplotlib.pyplot as plt

import numpy as np
from datetime import datetime, timedelta
import pickle
import netCDF4 as nc
import time

pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

# specify input
Ldir['gtagex'] = Ldir['gtag'] + '_lobio1'

if Ldir['parent'] == '/Users/PM5/Documents/':
    dt0 = datetime(2015,9,18)
    dt1 = datetime(2015,9,20)
    out_tag = 'test'
elif Ldir['parent'] == '/data1/parker/':
    dt0 = datetime(2015,1,1)
    dt1 = datetime(2015,12,31)
    out_tag = 'full'

# setup output location
out_dir0 = Ldir['parent'] + '/ptools_output/slow_slip/'
Lfun.make_dir(out_dir0, clean=False)
out_dir = out_dir0 + Ldir['gtagex'] + '_' + out_tag + '/'
Lfun.make_dir(out_dir, clean=True)

date_list = []
fn_list = []
dt = dt0
while dt <= dt1:
    date_list.append(dt.strftime('%Y.%m.%d'))
    dt = dt + timedelta(1)
for dl in date_list:
    f_string = 'f' + dl
    fn = (Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/'
        + f_string + '/low_passed.nc')
    fn_list.append(fn)

# calculate average bottom pressure
g = 9.8
tt = 0
NT = len(fn_list)
for fn in fn_list:
    ds = nc.Dataset(fn)
    zeta = ds['zeta'][:].squeeze()
    if tt ==0:
        G, S, T = zrfun.get_basic_info(fn)
        bp_arr = (0*zeta) * np.ones((NT,1,1))
        DA = G['DX'] * G['DY']
        DAm = np.ma.masked_where(zeta.mask, DA)
    rho = ds['rho'][:].squeeze() + 1000.
    # note that rho appears to be in situ density, not potential density 
    z_w = zrfun.get_z(G['h'], zeta, S, only_w=True)
    DZ = np.diff(z_w, axis=0)
    bp_arr[tt,:,:] = (g * rho * DZ).sum(axis=0) 
    tt += 1
bp_mean = np.mean(bp_arr, axis=0)

bp_anom = bp_arr - bp_mean

# plotting

plt.close('all')
cmap = plt.get_cmap(name='bwr')

tt = 0
vscale1 = 300
vscale2 = 300
for fn in fn_list:
    ds = nc.Dataset(fn)
    [T] = zrfun.get_basic_info(fn, getG=False, getS=False)
    fig = plt.figure(figsize=(12,8))
    
    bpa = bp_anom[tt, :, :]    
    # also make the anomaly away from the area mean at this time
    bpaa = bpa - ( (bpa*DAm).sum() / DAm.sum() )
    
    ax = fig.add_subplot(121)    
    cs = ax.pcolormesh(G['lon_psi'], G['lat_psi'], bpa[1:-1, 1:-1], vmin=-vscale1, vmax=vscale1, cmap='bwr')    
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    f_string = 'f' + datetime.strftime(T['tm'],'%Y.%m.%d')
    ax.set_title(f_string + ' Bottom $P^{\prime}$ [Pa]')
    
    ax = fig.add_subplot(122)    
    cs = ax.pcolormesh(G['lon_psi'], G['lat_psi'], bpaa[1:-1, 1:-1], vmin=-vscale2, vmax=vscale2, cmap='bwr')    
    fig.colorbar(cs)
    pfun.add_bathy_contours(ax, ds, txt=True)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    f_string = 'f' + datetime.strftime(T['tm'],'%Y.%m.%d')
    print('Working on ' + f_string)
    ax.set_title(' Bottom $P^{\prime\prime}$ [Pa]')
    
    save_string = ('000' + str(tt))[-3:]    
    plt.savefig(out_dir + 'p_prime_' + save_string + '.png')
    
    plt.close()
    
    tt += 1
