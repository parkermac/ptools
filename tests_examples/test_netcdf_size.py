"""
Code to test the size of a NetCDF file that has many masked or NaN values.

RESULT: it worked.  The full file was 7.1 MB and the masked one was 360 KB,
a factor of 20 smaller.  Note that the full file without compression is 8 MB,
consistent with 8 bytes for a single 64 bit float.

Without zlib=True there is no difference in size (both 8 MB).

With format='NETCDF3_64BIT_OFFSET' there is no difference in size (both 8 MB).
This is disappointing because that it what I use in ocn4.

From help(nc.Dataset):
**`data_model`**: `data_model` describes the netCDF
data model version, one of `NETCDF3_CLASSIC`, `NETCDF4`,
`NETCDF4_CLASSIC`, `NETCDF3_64BIT_OFFSET` or `NETCDF3_64BIT_DATA`.
Also, the default for "format" is NETCDF4.

In my current LO history files the format is NETCDF3 and the data_model
is NETCDF3_64BIT_OFFSET.  Also in the LO_ROMS makefile we do not appear
to ask it to use NETCDF4:


The results are basically identical when omitting complevel.
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
        
ncformat = 'NETCDF3_64BIT_OFFSET'

N = 1000 # grid size
n = 10 # edge band size
a = np.random.randn(N,N)

b = a.copy()
# This nan's out a large square, and keeps a band around the edges.
b[n:-n,n:-n] = np.nan

ds = nc.Dataset(out1_fn, 'w', data_model=ncformat)
ds.createDimension('x', N)
ds.createDimension('y', N)
#vv = ds.createVariable('a', float, ('x','y'), complevel=9, zlib=True)
vv = ds.createVariable('a', float, ('x','y'), zlib=True)
vv[:] = a
ds.close()

ds = nc.Dataset(out2_fn, 'w', data_model=ncformat)
ds.createDimension('x', N)
ds.createDimension('y', N)
#vv = ds.createVariable('a', float, ('x','y'), complevel=9, zlib=True)
vv = ds.createVariable('a', float, ('x','y'), zlib=True)
vv[:] = b
ds.close()

