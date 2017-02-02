#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 14:18:11 2017

@author: PM5

Code to plot the salinity variance field of Puget Sound
"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()
import zfun
import zrfun
pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt


which_home = os.environ.get("HOME") # This works even when called by cron.
if which_home == '/Users/PM5': # mac version
    in_dir = Ldir['parent'] + 'roms/output/salish_2006_4/'
elif which_home == '/home/parker': # fjord version
    in_dir = '/pmr3/pmraid1/daves/runs/salish_2006_4/OUT/'

ii = 4994
hh = ('000' + str(ii))[-4:]
fn = (in_dir + 'ocean_his_' + hh + '.nc')
G, S, T = zrfun.get_basic_info(fn)

x = G['lon_rho'][0,:]
y = G['lat_rho'][:,0]

x0 = -123.2
x1 = -122
y0 = 47
y1 = 48.4

i0 = zfun.find_nearest_ind(x, x0)
i1 = zfun.find_nearest_ind(x, x1)
j0 = zfun.find_nearest_ind(y, y0)
j1 = zfun.find_nearest_ind(y, y1)

h = G['h']

ds = nc.Dataset(fn)
eta = ds['zeta'][:].squeeze()

zr, zw = zrfun.get_z(h, eta, S)

DZ = np.diff(zw, axis = 0)
DZw = np.diff(zr, axis = 0)
DA = G['DX'] * G['DY']
DV = DZ * (DA * np.ones((S['N'], 1, 1)))

salt3d = ds['salt'][:].squeeze()
K = ds['AKs'][:].squeeze()

# vertically-averaged salinity     
sbar2d = np.sum(salt3d * DZ, axis = 0) / h
             
sbar3d = sbar2d * np.ones((S['N'], 1, 1))

sprime3d = salt3d - sbar3d

dsp_dz = np.diff(sprime3d, axis = 0) / DZw

mix2d = - 2  * np.sum(K[1:-1, :, :] * dsp_dz**2 * DZw, axis = 0)

# vertically-averaged LOCAL salinity variance
svar2d = np.sum(sprime3d**2 * DZ, axis = 0) / h

# volume-averaged salinity (over a limited volume)            
sbar1d =  (np.nansum(sbar3d[:, j0:j1, i0:i1] * DV[:, j0:j1, i0:i1]) /
                   np.nansum(DV[:, j0:j1, i0:i1]))
                   
# straining terms

xr = G['lon_rho']
yr = G['lat_rho']

u3d = ds['u'][:].squeeze()
ubar2d = ds['ubar'][:].squeeze()
ubar3d = ubar2d * np.ones((S['N'], 1, 1))
uprime3d = u3d - ubar3d

v3d = ds['v'][:].squeeze()
vbar2d = ds['vbar'][:].squeeze()
vbar3d = vbar2d * np.ones((S['N'], 1, 1))
vprime3d = v3d - vbar3d

sprime3d_u = (sprime3d[:, :, 1:] + sprime3d[:, :, :-1])/2
sprime3d_v = (sprime3d[:, 1:, :] + sprime3d[:, :-1, :])/2

DZ_u = (DZ[:, :, 1:] + DZ[:, :, :-1])/2
DZ_v = (DZ[:, 1:, :] + DZ[:, :-1, :])/2

DX = G['DX']
DY = G['DY']

DX_u = (DX[:, 1:] + DX[:, :-1])/2
DY_v = (DY[1:, :] + DY[:-1, :])/2
        
DX_u3 = DX_u  * np.ones((S['N'], 1, 1))
DY_v3 = DY_v * np.ones((S['N'], 1, 1))

hupspbar2d = np.sum(uprime3d * sprime3d_u * DZ_u, axis=0)
hvpspbar2d = np.sum(vprime3d * sprime3d_v * DZ_v, axis=0)

dsbardx2d = (sbar2d[:, 1:] - sbar2d[:, :-1]) / DX_u
dsbardy2d = (sbar2d[1:, :] - sbar2d[:-1, :]) / DY_v
             
dsprimedx3d = (sprime3d[:, :, 1:] - sprime3d[:, :, :-1]) / DX_u3
dsprimedy3d = (sprime3d[:, 1:, :] - sprime3d[:, :-1, :]) / DY_v3

strainx2d = 2 * hupspbar2d * dsbardx2d
strainy2d = 2 * hvpspbar2d * dsbardy2d

husp_triple = np.sum(uprime3d * sprime3d_u**2 * DZ_u, axis=0)
hvsp_triple = np.sum(vprime3d * sprime3d_v**2 * DZ_v, axis=0)

strain2d_triple = ( ((husp_triple[:, 1:] - husp_triple[:, :-1]) / DX[:, 1:-1])[1:-1, :]
                     + ((hvsp_triple[1:, :] - hvsp_triple[:-1, :]) / DY[1:-1, :])[:, 1:-1] )

#strainx2d_r = zfun.interp_scattered_on_plaid(xr.flatten(), yr.flatten(),
#    G['lon_u'][0, :], G['lat_u'][:, 0], strainx2d).reshape(xr.shape)
#strainy2d_r = zfun.interp_scattered_on_plaid(xr.flatten(), yr.flatten(),
#    G['lon_v'][0, :], G['lat_v'][:, 0], strainy2d).reshape(xr.shape)

strainx2d_r = (strainx2d[:, 1:] + strainx2d[:, :-1])/2
strainy2d_r = (strainy2d[1:, :] + strainy2d[:-1, :])/2

strain2d = strainx2d_r[1:-1, :] + strainy2d_r[:, 1:-1]

#strain2d = strain2d + strain2d_triple
                         
# plotting

plt.close('all')

fig = plt.figure(figsize=(18,10))

ax = fig.add_subplot(121)
vv = .005
cs = ax.pcolormesh(G['lon_psi'], G['lat_psi'],  strain2d,
      vmin=-vv, vmax = vv,  cmap='bwr')
ax.axis([x0, x1, y0, y1])
fig.colorbar(cs)
pfun.add_coast(ax)
pfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Vertically Integrated Straining (psu^2 m s-1)')

ax = fig.add_subplot(122)
cs = ax.pcolormesh(G['lon_psi'], G['lat_psi'],
    np.log10(-mix2d[1:-1, 1:-1] + 1e-10), vmin=-6, vmax = -2, cmap='rainbow')
ax.axis([x0, x1, y0, y1])
fig.colorbar(cs)
pfun.add_coast(ax)
pfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('log10[ Vertically Integrated (-) Mixing (psu^2 m s-1) ]')

plt.show()