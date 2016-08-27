# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 07:12:56 2016

@author: PM5

Code to experiment with adding a variable
to an existing NetCDF file.

"""

import netCDF4 as nc
import numpy as np

#%% check results
def check(fn):
    foo = nc.Dataset(fn)
    for vn in foo.variables:
        print(foo[vn])
        print(foo[vn][:])
        print('\n')
    foo.close()

fn1 = '/Users/PM5/Desktop/test1.nc'
fn2 = '/Users/PM5/Desktop/test2.nc'

ds1 = nc.Dataset(fn1, mode='w', format='NETCDF4')
# NETCDF3_CLASSIC also works,
ds1.createDimension('vec1', 4)
v_var = ds1.createVariable('Vec1', float, ('vec1',))
v_var[:] = [1,2,3,4]
ds1.close()
print('\nBefore')
check(fn1)

# open existing file.  Mode is 'a' for append.
ds2 = nc.Dataset(fn1, 'a')
# create new variable with dimensions of time (an existing dimension in the original file)
newvar = ds2.createVariable('newvar', np.float64, ('vec1',) )
# assign values
newvar[:] = np.ones( (len(ds2.dimensions['vec1']),) )
# close writes the file
ds2.close()

print('\nAfter')
check(fn1)

#%% This copies dimensions and variables - but I'm not sure when
# it would be useful.

if False:
    ds1 = nc.Dataset(fn1, mode='r')
    ds2 = nc.Dataset(fn2, mode='w')
    #Copy dimensions
    for dname, the_dim in ds1.dimensions.items():
        print(dname, len(the_dim))
        ds2.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
    # Copy variables
    for v_name, varin in ds1.variables.items():
        outVar = ds2.createVariable(v_name, varin.datatype, varin.dimensions)
        print(varin.datatype)
        # Copy variable attributes
        outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
        outVar[:] = varin[:]
    ds1.close()
    # add a variable to ds2
    ds2.createDimension('vec2', 3)
    v_var = ds2.createVariable('Vec2', float, ('vec2'))
    v_var[:] = [1,2,3]
    ds2.close()
    print('\nAfter')
    check(fn2)

