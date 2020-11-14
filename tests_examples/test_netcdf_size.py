"""
Code to test the size of a NetCDF file that has many masked or NaN values.

RESULT: it worked.  The full file was 7 MB and the masked one was 1.4 MB,
a factor of 5 smaller.  Note that the full file without compression is 8 MB,
consistent with 8 bytes for a single 64 bit float.  Also note that the masked version
only has 1% of the unmasked data, so the compression is nowhere near that.

Note: I'm not sure what complevel means.  Need to search around more.

"""

import netCDF4 as nc
import numpy as np

import os
out_dir = os.path.abspath('../../ptools_output/tests') + '/'
out1_fn = out_dir + 'f_nomask.nc'
out2_fn = out_dir + 'f_mask.nc'

# get rid of the old versions
for fn in [out1_fn, out2_fn]:
    try:
        os.remove(fn)
    except OSError:
        pass

N = 1000
n = 100
a = np.random.randn(N,N)

b = a.copy()
b[n:,n:] = np.nan

ds = nc.Dataset(out1_fn, 'w')
ds.createDimension('x', N)
ds.createDimension('y', N)
vv = ds.createVariable('a', float, ('x','y'), complevel=9, zlib=True)
vv[:] = a
ds.close()

ds = nc.Dataset(out2_fn, 'w')
ds.createDimension('x', N)
ds.createDimension('y', N)
vv = ds.createVariable('a', float, ('x','y'), complevel=9, zlib=True)
vv[:] = b
ds.close()

