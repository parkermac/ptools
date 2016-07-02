# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 15:34:17 2016

@author: PM5

Code to make a ROMS grid file.
"""

from importlib import reload
import gfun; reload(gfun)
G = gfun.gstart()

import numpy as np
import h5py
import netCDF4 as nc

import os

import Lfun
import zfun
import matfun


Lfun.make_dir(G['gdir'], clean=True)

fn = 'grid_m00_r00_s00_x00.nc'
out_fn = G['gdir'] + fn
print(50*'*')
print(out_fn)

# vectors to define the plaid grid
# start with cell corners (like an extended psi grid)

if G['gridname'] == 'test':
    # cascadia-like
    plon_vec = np.linspace(-127,-122,100)
    plat_vec = np.linspace(43,50,200)

elif G['gridname'] == 'ps':
    # Puget Sound
    plon_vec = np.linspace(-124,-122,400)
    plat_vec = np.linspace(47,49,600)

plon, plat = np.meshgrid(plon_vec, plat_vec)
ax_lims = (plon_vec[0], plon_vec[-1], plat_vec[0], plat_vec[-1])

# specify topography files to use
t_dir = G['dir0'] + 'tools_data/geo_data/topo/'
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
z = np.nan * lon

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

    z_flat = zfun.interp_scattered_on_plaid(lon.flatten(), lat.flatten(),
                                         tlon_vec, tlat_vec, tz)
    z_part = z_flat.reshape((NR, NC))

    # put good values of z_part in z
    z[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]

#%% Adjust zero of the bathymetry to account for the fact that mean sea level
# is somewhat higher than NAVD88.
z = z - 1.06

#%% save the output to NetCDF

# get rid of old version
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist

# create new NetCDF file
foo = nc.Dataset(out_fn, 'w')

# create dimensions
M, L = lon.shape # use ROMS teminology
size_dict = {'rho': (M, L),
             'u': (M, L-1),
             'v': (M-1, L),
             'psi': (M-1, L-1),
             'psi_ex': (M+1, L+1)}
tag_list = ['rho', 'u', 'v', 'psi', 'psi_ex']
for tag in tag_list:
    foo.createDimension('eta_'+tag, size_dict[tag][0])
    foo.createDimension('xi_'+tag, size_dict[tag][1])

# create variables
lon_var = dict()
lat_var = dict()
for tag in tag_list:
    lat_var[tag] = foo.createVariable('lat_'+tag, float, ('eta_'+tag, 'xi_'+tag))
    lon_var[tag] = foo.createVariable('lon_'+tag, float, ('eta_'+tag, 'xi_'+tag))
z_var = foo.createVariable('z', float, ('eta_rho', 'xi_rho'))
mask_var = foo.createVariable('mask_rho', int, ('eta_rho', 'xi_rho'))
dx_var = foo.createVariable('dx', float, ('eta_rho', 'xi_rho'))
dy_var = foo.createVariable('dy', float, ('eta_rho', 'xi_rho'))
pm_var = foo.createVariable('pm', float, ('eta_rho', 'xi_rho'))
pn_var = foo.createVariable('pn', float, ('eta_rho', 'xi_rho'))

# create other grids
lon_dict = dict()
lat_dict = dict()

lon_dict['rho'] = lon
lat_dict['rho'] = lat

lon_dict['psi_ex'] = plon
lat_dict['psi_ex'] = plat

lon_dict['u'] = lon[:, :-1] + np.diff(lon, axis=1)/2
lat_dict['u'] = lat[:, :-1]

lon_dict['v'] = lon[:-1, :]
lat_dict['v'] = lat[:-1, :] + np.diff(lat, axis=0)/2

lon_dict['psi'] = plon[1:-1, 1:-1]
lat_dict['psi'] = plat[1:-1, 1:-1]

# create dx, dy and pm, pn
ulon = plon[:-1, :]
vlat = plat[:, :-1]
R = zfun.earth_rad(np.mean(plat[:,1]))
dx = R * np.cos(np.pi*lat/180) * (np.pi*np.diff(ulon, axis=1)/180)
dy = R * (np.pi*np.diff(vlat, axis=0)/180)

# add data to fields
for tag in tag_list:
    lon_var[tag][:] = lon_dict[tag]
    lat_var[tag][:] = lat_dict[tag]

z_var[:] = z
mask_rho = np.ones((M, L)) # start with all ones (unmasked for ROMS)
mask_var[:] = mask_rho
dx_var[:] = dx
dy_var[:] = dy

foo.close()
