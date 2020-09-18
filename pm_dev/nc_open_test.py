"""
Code to test faster ways to open NetCDF files.

inspired by:
https://github.com/Unidata/netcdf-c/issues/489

"""

import netCDF4 as nc
import xarray as xr
from time import time
import numpy as np

indir = '/Users/pm8/Documents/LiveOcean_roms/output/cas6_v3_lo8b/f2019.07.04/'

tt0 = time()

NT = 25
for ii in range(1,NT+1):
    fn = indir + 'ocean_his_' + ('0000' + str(ii))[-4:] + '.nc'
    ds = nc.Dataset(fn)
    s = ds['salt'][0,:,:,:]
    if ii == 1:
        N, M, L = s.shape
        Maskr = ds['mask_rho'][:] == 1 # True over water
        Maskr3 = np.tile(Maskr.reshape(1, M, L),[N,1,1])
    ss = s[Maskr3].data
    ds.close()
    
print('Took %0.1f sec for %d opens with netCDF4' % (time()-tt0, NT))

tt0 = time()

NT = 25
for ii in range(1,NT+1):
    if ii == 1:
        mm = Maskr3.data
    fn = indir + 'ocean_his_' + ('0000' + str(ii))[-4:] + '.nc'
    A = xr.open_dataset(fn)
    a = A.salt.squeeze().data
    aa = a[mm]
    A.close()
    
print('Took %0.1f sec for %d opens with xarray' % (time()-tt0, NT))