"""
Create a tidally averaged surface Bernoulli function
from a LiveOcean layer extraction.
"""

# setup
import netCDF4 as nc
import numpy as np
from time import time

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
import zfun

import subprocess

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

Ldir = Lfun.Lstart()
in_fn = Ldir['LOo'] + 'layer/cas6_v3_lo8b_2019.06.01_2019.08.31/surface_hourly.nc'
# a total of 92 days, so we can form 90 tidal averages I believe


# make a place to save output
outdir = os.path.abspath('../../ptools_output/bernoulli_average') + '/'
Lfun.make_dir(outdir)

out_fn = outdir+'surface_lp.nc'

# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist

cmd_list = ['ncks', '-v', 'lon_psi,lat_psi,lon_rho,lat_rho,ocean_time,zeta', '-d','ocean_time,0,0', '-O', in_fn, out_fn]
# NOTE this should work with '-d','ocean_time,0' because this should retain the time dimension

tt0 = time()
ret1 = subprocess.call(cmd_list)
print('Time initialize using ncks = %0.1f sec' % (time() - tt0))

ds_out = nc.Dataset(out_fn, 'a')
ds_out.createVariable('uuvv', float, ('ocean_time', 'eta_rho', 'xi_rho'))
ds_out.createVariable('u', float, ('ocean_time', 'eta_rho', 'xi_rho'))
ds_out.createVariable('v', float, ('ocean_time', 'eta_rho', 'xi_rho'))

f = zfun.godin_shape()
NT = len(f)

ds = nc.Dataset(in_fn)
ot = ds['ocean_time'][:]
nt = len(ot)

day_list = [item for item in range(90)]
# valid values are 0-89 so range(90) works
tt0 = time()
for day in day_list:
    print(day)
    i0 = 1 + day*24
    ii = 0
    for tt in range(i0,i0+NT):
        if tt == i0:
            t = ot[tt]/NT
            zr = ds['zeta'][tt,:,:].squeeze()*f[ii]
            ur = ds['u'][tt,:,:].squeeze()
            vr = ds['v'][tt,:,:].squeeze()
            u = ur*f[ii]
            v = vr*f[ii]
            uuvv = (ur*ur + vr*vr)*f[ii]
            uuvv[uuvv.mask] = 0
            uuvv[zr.mask] = np.nan
        else:
            t = t + ot[tt]/NT
            zr = zr + ds['zeta'][tt,:,:].squeeze()*f[ii]
            ur = ds['u'][tt,:,:].squeeze()
            vr = ds['v'][tt,:,:].squeeze()
            u = u + ur*f[ii]
            v = v + vr*f[ii]
            uuvv1 = (ur*ur + vr*vr)*f[ii]
            uuvv1[uuvv1.mask] = 0
            uuvv = uuvv + uuvv1
        ii += 1
    ds_out['zeta'][day,:,:] = zr
    ds_out['u'][day,:,:] = u
    ds_out['v'][day,:,:] = v
    ds_out['uuvv'][day,:,:] = uuvv
    ds_out['ocean_time'][day] = t
    
print('Time to add averages = %0.1f sec' % (time() - tt0))

for vn in ds_out.variables:
    print('')
    print(ds_out[vn])
    
ds_out.close()
ds.close()

