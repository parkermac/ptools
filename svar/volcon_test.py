"""
Code to test calculation of volume conservation.
"""

import os
import sys
import netCDF4 as nc
import numpy as np

pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import zfun
import zrfun

dir00 = '/Users/pm7/Documents/LiveOcean_roms/output/aestus1_A1_ae1/'
dir0 = dir00 + 'f2013.03.02/'

fn_h0 = dir0 + 'ocean_his_0001.nc'
fn_h1 = dir0 + 'ocean_his_0002.nc'
fn_a = dir0 + 'ocean_avg_0001.nc'

ds_h0 = nc.Dataset(fn_h0)
ds_h1 = nc.Dataset(fn_h1)
ds_a = nc.Dataset(fn_a)

G, S, T0 = zrfun.get_basic_info(fn_h0)
T1 = zrfun.get_basic_info(fn_h1, only_T=True)
DT = T1['ocean_time'] - T0['ocean_time']
DA = G['DX'][1:-1, 1:-1] * G['DY'][1:-1, 1:-1]

zeta0 = ds_h0['zeta'][:].squeeze()
zeta1 = ds_h1['zeta'][:].squeeze()

Huon = ds_a['Huon'][:].squeeze()
Hvom = ds_a['Hvom'][:].squeeze()

fx = Huon.sum(axis=0)
fy = Hvom.sum(axis=0)

transport_in = -(np.diff(fx, axis=1)[1:-1, :]
    + np.diff(fy, axis=0)[:, 1:-1])

dvol_dt = DA *(zeta1[1:-1, 1:-1] - zeta0[1:-1, 1:-1])/DT

err = dvol_dt - transport_in

rel_err = err / dvol_dt

ds_h0.close()
ds_h1.close()
ds_a.close()

import matplotlib.pyplot as plt
plt.close('all')
fig = plt.figure(figsize=(12,7))
ax = fig.add_subplot(111)
remax = .05
cs = ax.pcolormesh(G['lon_psi'], G['lat_psi'], rel_err, vmin=-remax, vmax = remax, cmap='jet')
fig.colorbar(cs)
plt.show()

