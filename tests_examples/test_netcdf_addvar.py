# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 07:12:56 2016

@author: PM5

Experiment with adding a variable to an existing NetCDF file.

"""

import netCDF4 as nc
import numpy as np
import os

dir0 = os.environ.get('HOME') + '/Desktop/'

# function to check results
def check(fn, hstr):
    print('\n** ' + hstr + '\n' + 60*'=')
    ds = nc.Dataset(fn)
    for vn in ds.variables:
        #print(ds[vn])
        print(vn + str(ds[vn].dimensions) + ': ' + str(ds[vn][:]))
    ds.close()

# output file name
fn = dir0 + 'run test_n test1.nc'
try:
    os.remove(fn)
except FileNotFoundError:
    pass

# size to use for data
N = 4

hstr = 'Create the original file.'
ds = nc.Dataset(fn, 'w', format='NETCDF4')
#ds = nc.Dataset(fn, 'w', format='NETCDF3_CLASSIC')
# NETCDF3_CLASSIC also works
ds.createDimension('d1', None)
vv = ds.createVariable('data1', float, 'd1')
vv[:] = range(N)
ds.close()
check(fn, hstr)

hstr = 'Add a new variable using an existing dimension.'
ds = nc.Dataset(fn, 'a')
vv = ds.createVariable('data2', float, 'd1')
NN = len(ds.dimensions['d1'])
vv[:] = np.ones(NN)
ds.close()
check(fn, hstr)

hstr = 'Append data to the end of an existing variable.'
ds = nc.Dataset(fn, 'a')
ds.variables['data2'][N: 2*N] = range(N)
ds.close()
check(fn, hstr)

hstr = 'Add a dimension and use it to create a new variable.'
ds = nc.Dataset(fn, 'a')
ds.createDimension('d3', N)
vv = ds.createVariable('data3', float, 'd3')
vv[:] = range(N)
ds.close()
check(fn, hstr)

hstr = 'Replace a variable with zeros.\n  NOTE that you cannot delete a variable.'
ds = nc.Dataset(fn, 'a')
ds['data3'][:] = 0.
ds.close()
check(fn, hstr)


