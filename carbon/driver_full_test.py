"""
Code to test using CO2SYS.m

Test of doing carbon calc for one full history file.

RESULT:

TIME to gather and process input fields:
 -- 21.6 seconds

 TIME to do calculation:
 -- run CO2SYS 313.9 s

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()
import zrfun

import carbon_fun as cfun
from importlib import reload
reload(cfun)

import scipy.io as io
import netCDF4 as nc
import subprocess
import numpy as np
import time
import seawater as sw

def fillit(a):
    # ensures a is an array with nan's for masked values
    # instead of a masked array
    if isinstance(a, np.ma.MaskedArray):
        a = a.filled(np.nan)
    return a

# (1) Get and package input fields
# file to work on
fn = (Ldir['roms'] + '/output/' + 'cas4_v2_lo6biom/f2018.09.29/ocean_his_0001.nc')

tt0 = time.time()
ds = nc.Dataset(fn)
G, S, T = zrfun.get_basic_info(fn)


# extract needed info from history file
v_dict = dict()
for vn in ['alkalinity', 'TIC', 'salt', 'temp','rho']:
    v = ds[vn][:]
    v = fillit(v)
    v_dict[vn] = v
            
# create depth, pressure, and in situ temperature
h = ds['h'][:]
h = fillit(h)
lat = G['lat_rho'][:]
z_rho = zrfun.get_z(h, 0*h, S, only_rho=True)
depth = -z_rho
pres = sw.pres(depth, lat)
v_dict['pres'] = pres
temp = sw.ptmp(v_dict['salt'], v_dict['temp'], 0, v_dict['pres'])
v_dict['temp'] = temp

# convert from umol/L to umol/kg using in situ dentity
v_dict['alkalinity'] = 1000 * v_dict['alkalinity'] / (v_dict['rho'] + 1000)
v_dict['TIC'] = 1000 * v_dict['TIC'] / (v_dict['rho'] + 1000)

# testing
# v_dict['alkalinity'][v_dict['alkalinity']<100] = 100
# v_dict['TIC'][v_dict['TIC']<100] = 100

# clean up
v_dict.pop('rho') # no longer needed, so don't pass to worker
ds.close()

print('\nTIME to gather and process input fields:')
print(' -- %0.1f seconds' % (time.time() - tt0))
    
# (2) run CO2SYS
tempdir = '../../ptools_output/carbon/temp/'
PH, OM = cfun.get_carbon(v_dict, tempdir, print_info=True)

