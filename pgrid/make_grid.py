# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 15:34:17 2016

@author: PM5

Code to initialize the creation of a ROMS grid file.

NOTE: the gridname is set in gfun.gstart().

Throughout this code I try to use ROMS naming conventions, except that
when manipulating or plotting I refer to [lon,lat]_rho as [lon,lat],
and [lon,lat]_psi_ex as [plon,plat].
"""

from importlib import reload
import gfun; reload(gfun)
G = gfun.gstart()

import numpy as np
import h5py
import netCDF4 as nc
import seawater as sw

import os

import Lfun
import zfun
import matfun

#%%
Lfun.make_dir(G['gdir'], clean=True)

fn = 'grid_m00_r00_s00_x00.nc'
out_fn = G['gdir'] + fn
print(50*'*')
print(out_fn)

# vectors to define the plaid grid
# start with cell corners (like an extended psi grid)

def simple_grid(aa, res):
    dlat = aa[3] - aa[2]
    dlon = aa[1] - aa[0]
    mean_lat = np.mean(aa[2:])
    earth_rad = zfun.earth_rad(mean_lat)
    dlatm = earth_rad * np.pi * dlat/180
    dlonm = earth_rad * np.cos(np.pi*mean_lat/180) * np.pi * dlon/180
    nx = int(np.ceil(dlonm/res))
    ny = int(np.ceil(dlatm/res))
    plon_vec = np.linspace(aa[0], aa[1], nx)
    plat_vec = np.linspace(aa[2], aa[3], ny)
    print('<< ' + G['gridname'] + ' >>')
    print('   nx = %d' % (nx))
    print('   ny = %d' % (ny))
    return plon_vec, plat_vec

def stretched_grid(lon_list, x_res_list, lat_list, y_res_list):
    """
    Input:
    The four lists (units = degrees and meters) are points that define
    segmented vectors along which the resolution changes linearly.
    """
    plon_list = []
    plat_list = []
    if (len(lon_list) != len(x_res_list)) or (len(lat_list) != len(y_res_list)):
        print('Lists must be the same length')
        return np.array(plon_list), np.array(plat_list)

    lon_vec = np.array(lon_list)
    x_res_vec = np.array(x_res_list)
    lat_vec = np.array(lat_list)
    y_res_vec = np.array(y_res_list)

    R = zfun.earth_rad(np.mean(lat_vec))
    clat = np.cos(np.pi*np.mean(lat_vec)/180)

    lon = lon_list[0]
    plon_list.append(lon)
    while lon <= lon_list[-1]:
        i0, i1, fr = zfun.get_interpolant(np.array([lon]), lon_vec)
        xres = (1-fr)*x_res_vec[i0] + fr*x_res_vec[i1]
        dlon = 180 * xres / (clat * R * np.pi)
        lon = lon + dlon
        plon_list.append(lon)
    lat = lat_list[0]
    plat_list.append(lat)
    while lat <= lat_list[-1]:
        i0, i1, fr = zfun.get_interpolant(np.array([lat]), lat_vec)
        yres = (1-fr)*y_res_vec[i0] + fr*y_res_vec[i1]
        dlat = 180 * yres / (R * np.pi)
        lat = lat + dlat
        plat_list.append(lat)
    return np.array(plon_list), np.array(plat_list)

if G['gridname'] == 'cascadia2':
    # cascadia-like
    aa = [-127.4, -122, 43, 50]
    res = 5000 # target resolution (m)
    plon_vec, plat_vec = simple_grid(aa, res)
    
if G['gridname'] == 'salish1':
    # focus on Puget Sound
    aa = [-124, -122, 46.8, 49]
    res = 300 # target resolution (m)
    plon_vec, plat_vec = simple_grid(aa, res)
    
elif G['gridname'] == 'bigSalish1':
    aa = [-130, -122, 43, 51.5]
    res = 2000 # target resolution (m)
    plon_vec, plat_vec = simple_grid(aa, res)
    
elif G['gridname'] == 'bigSalish2':
    lon_list = [-130, -127, -122]
    x_res_list = [5000, 500, 500]
    lat_list = [42, 44, 47, 51, 51.5]
    y_res_list = [5000, 2000, 500, 500, 2000]
    plon_vec, plat_vec = stretched_grid(lon_list, x_res_list,
                                        lat_list, y_res_list)
elif G['gridname'] == 'aestus1':
    lon_list = [-1, 0, 1, 2, 3]
    x_res_list = [5000, 1000, 1000, 5000, 5000]
    lat_list = [44, 44.9, 45.1, 46]
    y_res_list = [5000, 1000, 1000, 5000]
    plon_vec, plat_vec = stretched_grid(lon_list, x_res_list,
                                        lat_list, y_res_list)

plon, plat = np.meshgrid(plon_vec, plat_vec)
ax_lims = (plon_vec[0], plon_vec[-1], plat_vec[0], plat_vec[-1])

# specify topography files to use
t_dir = G['dir0'] + 'tools_data/geo_data/topo/'
# list of topo files: coarsest to finest
#t_list = ['smith_sandwell/pnw_smithsand.mat',
#          'cascadia/cascadia_gridded.mat',
#         'psdem/PS_183m.mat']
t_list = ['srtm15/topo15.grd',
          'cascadia/cascadia_gridded.mat',
         'psdem/PS_183m.mat',
         'ttp_patch/TTP_Regional_27m_patch.mat']

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
    
def load_bathy2(t_fn, lon_vec, lat_vec):
    # load a section of the new NetCDF Smith-Sandwell bathy
    ds = nc.Dataset(t_fn)
    Lon = ds['lon'][:]
    Lat = ds['lat'][:]
    i0 = zfun.find_nearest_ind(Lon, lon_vec[0])
    i1 = zfun.find_nearest_ind(Lon, lon_vec[-1])
    j0 = zfun.find_nearest_ind(Lat, lat_vec[0])
    j1 = zfun.find_nearest_ind(Lat, lat_vec[-1])
    tlon_vec = Lon[i0-2:i1+3]
    tlat_vec = Lat[j0-2:j1+3]
    tz = ds['z'][j0-2:j1+3, i0-2:i1+3]
    ds.close()
    return tlon_vec, tlat_vec, tz

if G['gridname'] == 'aestus1':
    # make grid and bathymetry by hand
    z = np.zeros(lon.shape)
    x, y = zfun.ll2xy(lon, lat, 0, 45)
    zshelf = x * 1e-3
    zestuary = -20 + 20*x/1e5 + 20/(1e4)*np.abs(y)
    z = zshelf
    mask = zestuary < z
    z[mask] = zestuary[mask]

else:
    # add bathymetry automatically from files
    for t_file in t_list:
        t_fn = t_dir + t_file
        print('\nOPENING BATHY FILE: ' + t_file)
        if t_file == 'srtm15/topo15.grd':
            tlon_vec, tlat_vec, tz = load_bathy2(t_fn, lon_vec, lat_vec)
        else:
            tlon_vec, tlat_vec, tz = load_bathy(t_fn)
        z_flat = zfun.interp_scattered_on_plaid(lon.flatten(), lat.flatten(),
                                             tlon_vec, tlat_vec, tz)
        z_part = z_flat.reshape((NR, NC))
        # put good values of z_part in z
        z[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]
    # Adjust zero of the bathymetry to account for the fact that mean sea level
    # is somewhat higher than NAVD88.
    z = z - 1.06

#%% save the output to NetCDF

# get rid of old version
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist

# create new NetCDF file
foo = nc.Dataset(out_fn, 'w', format='NETCDF3_CLASSIC')

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
mask_var = dict()
for tag in tag_list:
    lat_var[tag] = foo.createVariable('lat_'+tag, float, ('eta_'+tag, 'xi_'+tag))
    lon_var[tag] = foo.createVariable('lon_'+tag, float, ('eta_'+tag, 'xi_'+tag))
    mask_var[tag] = foo.createVariable('mask_'+tag, float, ('eta_'+tag, 'xi_'+tag))
h_var = foo.createVariable('h', float, ('eta_rho', 'xi_rho'))
f_var = foo.createVariable('f', float, ('eta_rho', 'xi_rho'))
pm_var = foo.createVariable('pm', float, ('eta_rho', 'xi_rho'))
pn_var = foo.createVariable('pn', float, ('eta_rho', 'xi_rho'))
xl_var = foo.createVariable('xl', float)
el_var = foo.createVariable('el', float)
sph_var = foo.createVariable('spherical', 'c')

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
    # start with all ones (unmasked for ROMS)
    mask_var[tag][:] = np.ones_like(lon_dict[tag], dtype=int)

h_var[:] = -z
pm_var[:] = 1/dx
pn_var[:] = 1/dy
f_var[:] = sw.f(lat_dict['rho'])
xl_var[:] = dx[0,:].sum()
el_var[:] = dy[:,0].sum()
sph_var[:] = 'T'

foo.close()
