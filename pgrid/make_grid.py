# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 15:34:17 2016

@author: PM5

Code to make a ROMS grid file.
"""

from importlib import reload
import gfun; reload(gfun)
gridname, dir0, pgdir, gdir = gfun.gstart()

import numpy as np
import h5py
import netCDF4 as nc

import os

import Lfun
import zfun
import matfun

Lfun.make_dir(gdir, clean=True)

fn = 'grid_m00_r00_s00_x00.nc'
out_fn = gdir + fn
print(50*'*')
print(out_fn)

# vectors to define the plaid grid
# start with cell corners (like an extended psi grid)

if gridname == 'test':
    # cascadia-like
    plon_vec = np.linspace(-127,-122,100)
    plat_vec = np.linspace(43,50,200)

elif gridname == 'ps':
    # Puget Sound
    plon_vec = np.linspace(-124,-122,400)
    plat_vec = np.linspace(47,49,600)

plon, plat = np.meshgrid(plon_vec, plat_vec)
ax_lims = (plon_vec[0], plon_vec[-1], plat_vec[0], plat_vec[-1])

# specify topography files to use
t_dir = dir0 + 'tools_data/geo_data/topo/'
# list of topo files: coarsest to finest
t_list = ['smith_sandwell/pnw_smithsand.mat',
          'cascadia/cascadia_gridded.mat',
         'psdem/PS_183m.mat']

# then make box centers
lon_vec = plon_vec[:-1] + np.diff(plon_vec)/2
lat_vec = plat_vec[:-1] + np.diff(plat_vec)/2
lon, lat = np.meshgrid(lon_vec, lat_vec)
NR, NC = lon.shape

# initialize the final bathymetry array
TZ = np.nan * lon

# coastline
c_dir = dir0 + 'tools_data/geo_data/coast/'
c_file = 'pnw_coast_combined.mat'
c_fn = c_dir + c_file
cmat = matfun.loadmat(c_fn)

def how_many_nans(arr):
    print(' Array: npts=' + str(arr.size) + ' nnan=' +
          str(sum(np.isnan(arr).flatten())))

def load_bathy(t_fn):
    try:
        a = h5py.File(t_fn)
        aa = dict()
        print(' using h5py')
        for item in a.keys():
            aa[item] = a[item][:]
            #print('    ' + item + ' ' + str(aa[item].shape))
        # reshaping because of how things are packed
        tlon_vec = a['lon'][:].flatten()
        tlat_vec = a['lat'][:].flatten()
        tz = a['z'][:].T
        #tmask = a['mask'][:].T
        # The field "mask" was created in the topo processing
        # as a matrix of ints with 1=ocean, 0=land.
        # We will supersede it here with our own masking.
        a.close()
    except:
        # for some reason h5py does not work with, for example,
        # psdem/PS_183m.mat
        tmat = matfun.loadmat(t_fn)
        print(' using matfun')
        # reshaping because of how things are packed
        tlon_vec = tmat['lon'][:].flatten()
        tlat_vec = tmat['lat'][:].flatten()
        tz = tmat['z'][:]
    return tlon_vec, tlat_vec, tz

for t_file in t_list:
    t_fn = t_dir + t_file
    print('\nOPENING BATHY FILE: ' + t_file)
    tlon_vec, tlat_vec, tz = load_bathy(t_fn)

    TZv = zfun.interp_scattered_on_plaid(lon.flatten(), lat.flatten(),
                                         tlon_vec, tlat_vec, tz)
    TZpart = TZv.reshape((NR, NC))

    # put good values of TZpart in TZ
    TZ[~np.isnan(TZpart)] = TZpart[~np.isnan(TZpart)]


#%% save the output to NetCDF

# get rid of old version
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist

# create new NetCDF file
M, L = lon.shape # use ROMS teminology
foo = nc.Dataset(out_fn, 'w')
foo.createDimension('eta_rho', M)
foo.createDimension('xi_rho', L)
foo.createDimension('y1', M + 1)
foo.createDimension('x1', L + 1)
lat_var = foo.createVariable('lat_rho', float, ('eta_rho', 'xi_rho'))
lon_var = foo.createVariable('lon_rho', float, ('eta_rho', 'xi_rho'))
# we call these "_psi_ex" becasue they extend one point beyond
# that of the standard ROMS psi grid
plat_var = foo.createVariable('lat_psi_ex', float, ('y1', 'x1'))
plon_var = foo.createVariable('lon_psi_ex', float, ('y1', 'x1'))
z_var = foo.createVariable('z', float, ('eta_rho', 'xi_rho'))
mask_var = foo.createVariable('mask_rho', int, ('eta_rho', 'xi_rho'))
lat_var[:] = lat
lon_var[:] = lon
plat_var[:] = plat
plon_var[:] = plon
z_var[:] = TZ
mask_rho = np.ones((M, L)) # start with all ones (unmasked for ROMS)
mask_var[:] = mask_rho
foo.close()
