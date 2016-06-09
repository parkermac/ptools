# -*- coding: utf-8 -*-
"""
Code to test nco_class.

"""
from importlib import reload
import nco_class
reload(nco_class)

dir0 = '/Users/PM5/Documents/LiveOcean_output/cascadia1_base/'
fn = dir0 + 'f2015.11.12/ocn/Data/t3d.nc'

#dir0 = '/Users/PM5/Documents/LiveOcean_output/cascadia1_base/'
#fn = dir0 + 'f2016.01.29/ocn/ocean_bry.nc'

#dir0 = '/Users/PM5/Documents/LiveOcean_roms/output/cascadia1_base_lo1/'
#fn = dir0 + 'f2016.01.29/ocean_his_0005.nc'

nco = nco_class.Nco(fn)

print('\n' + nco.fn)

nco.vlist()

#s = nco.ds['salt_west'][:]
#print('\nsalt_west')
#print(+ s[:,-1,100])
#
#t = nco.ds['temp_west'][:]
#print('\nsalt_west')
#print(t[:,-1,100])

a = nco.ds['t3d'][:]

print(a[:,-1,60,0])