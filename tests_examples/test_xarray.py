"""
Test of combining history files using xarray.

Typical commands:

z = a.zeta.mean(axis=0).values
s = a.salt.mean(axis=0).values
# interesting fact: it converts to correct datetimes!
# look at a.ocean_time.values

RESULT (7 3-D variables at one x,y for 49 times (2 days)):
extract using xarray = 0.26 sec
extract using netCDF4 = 0.17 sec
(on my mac)
** no apparent advantage using xarray in this test,
although there may be other times where lazy computation is desireable
or maybe when we use a full year of history files...?

This seems so fast - why does a mooring extraction take many hours for a year
on boiler?  The test above suggest it should take like a minute!
"""

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from time import time
import netCDF4 as nc

# this addresses some xarray plot warning
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

import sys, os
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart()

# get list of history files
Ldir['gtagex'] = 'cas6_v3_lo8b'
fn_list = Lfun.get_fn_list('hourly', Ldir, '2019.07.04', '2019.07.05')

# benchmark two ways of doing a mooring extraction
vn_list = ['u','v','salt','temp','oxygen','NO3','phytoplankton']

a = xr.open_mfdataset(fn_list, combine='nested', concat_dim='ocean_time',
        data_vars='minimal', coords='minimal', compat='override')
        # need combine to work with concat_dim
tt0 = time()
xr_dict = {}
for vn in vn_list:
     xr_dict[vn]= a[vn][:,:,10,10].values
print('extract using xarray = %0.2f sec' % (time()-tt0))

tt0 = time()
nc0 = np.nan * np.ones(xr_dict['u'].shape)
nc_dict = {}
for vn in vn_list:
    nc_dict[vn] = nc0.copy()
    
counter = 0
for fn in fn_list:
    ds = nc.Dataset(fn)
    for vn in vn_list:
        nc_dict[vn][counter,:] = ds[vn][0,:,10,10]
    ds.close()
    counter += 1
print('extract using netCDF4 = %0.2f sec' % (time()-tt0))

