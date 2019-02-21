"""
Code to test using CO2SYS.m

- Read in a roms history file.
- Save selected fields as matlab .mat files.
- Run a matlab subprocess to process these files.
- Read the results back into python and plot.

RESULT: CO2SYS takes 14 sec to run for one layer of the cas4 grid, on my mac.
All other tasks take minimul time.
"""

import scipy.io as io
import netCDF4 as nc
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import time
import seawater as sw

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()
import zrfun

fn = ('/Users/pm7/Documents/LiveOcean_roms/output/' +
    'cas4_v2_lo6biom/f2018.09.29/ocean_his_0001.nc')
ds = nc.Dataset(fn)

tt0 = time.time()
# extract needed info from history file
v_dict = dict()
for vn in ['alkalinity', 'TIC', 'salt', 'temp']:
    v = ds[vn][0,0,:,:].squeeze()
    v_dict[vn] = v
    vv = v.data
    vv[v.mask] = np.nan
    
    name = ds[vn].long_name
    try:
        units = ds[vn].units
    except AttributeError:
        units = ''
    vmax = np.nanmax(vv)
    vmin = np.nanmin(vv)
    v_dict[vn] = vv
    print('%25s (%25s) max = %6.1f min = %6.1f' %
        (name, units, vmax, vmin))
print('')
print('gather data %0.1f s' % (time.time()-tt0))
sys.stdout.flush()

tt0 = time.time()
# create potential density (note temp = potential temp)
pden = sw.dens(v_dict['salt'], v_dict['temp'], 0)
# convert from umol/L to umol/kg
v_dict['alkalinity'] = 1000 * v_dict['alkalinity'] / pden
v_dict['TIC'] = 1000 * v_dict['TIC'] / pden
print('pre-process data %0.1f s' % (time.time()-tt0))

tt0 = time.time()        
# save to temporary file for matlab
dir0 = '../../ptools_output/carbon/temp/'
Lfun.make_dir(dir0, clean=True)
in_fn = dir0 + 'input.mat'
io.savemat(in_fn,v_dict)
print('save data to .mat %0.1f s' % (time.time()-tt0))
sys.stdout.flush()

tt0 = time.time()
# run the matlab code
run_cmd = [Ldir['which_matlab'], "-nodisplay", "-r", "worker()", "&"]
proc = subprocess.run(run_cmd, stdout=subprocess.PIPE)
print('run CO2SYS %0.1f s' % (time.time()-tt0))
sys.stdout.flush()

tt0 = time.time()
# load the output of CO2SYS
PH = io.loadmat(dir0 + 'PH.mat')['PH']
OM = io.loadmat(dir0 + 'OM.mat')['OM']
print('load mat files to python %0.1f s' % (time.time()-tt0))
sys.stdout.flush()

plt.close('all')
fig = plt.figure(figsize=(16,8))

ax = fig.add_subplot(121)
cs = plt.pcolormesh(PH, vmin=7, vmax=8.5, cmap='Spectral')
plt.colorbar(cs)
ax.set_title('PH')

ax = fig.add_subplot(122)
cs = plt.pcolormesh(OM, vmin=0, vmax=3, cmap='coolwarm_r')
plt.colorbar(cs)
ax.set_title('Omega')

plt.show()
