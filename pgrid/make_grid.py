# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 15:34:17 2016

@author: PM5

Code to make a ROMS grid file.
"""

dir0 = '/Users/PM5/Documents/'

import numpy as np
import h5py
#import scipy.interpolate as intp
import numpy.ma as ma

#from warnings import filterwarnings
#filterwarnings('ignore') # skip some warning messages

import os
import sys
alp = os.path.abspath(dir0 +'LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import zfun
import matfun

outdir = dir0 +'ptools_output/pgrid/'

# name for the grid
gname_base = 'hires'
gname_tag = '_m00_r00_s00_x00'
gname = gname_base + gname_tag + '.nc'
out_fn = outdir + gname
print(50*'*')
print(out_fn)

# vectors to define the plaid grid
# start with cell corners (like an extended psi grid)
plon_vec = np.linspace(-127,-122,300)
plat_vec = np.linspace(43,50,600)
plon, plat = np.meshgrid(plon_vec, plat_vec)
ax_lims = (plon_vec[0], plon_vec[-1], plat_vec[0], plat_vec[-1])

# specify topography files to use
t_dir = dir0 + 'tools_data/geo_data/topo/'
# list files coarsest to finest
t_list = ['smith_sandwell/pnw_smithsand.mat',
          'cascadia/cascadia_gridded.mat',
         'psdem/PS_183m.mat']

# then make box centers
lon_vec = plon_vec[:-1] + np.diff(plon_vec)/2
lat_vec = plat_vec[:-1] + np.diff(plat_vec)/2
lon, lat = np.meshgrid(lon_vec, lat_vec)

# initialize the final bathymetry array
TZ0 = np.nan * lon

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

    NR, NC = lon.shape
    TZv = zfun.interp_scattered_on_plaid(lon.flatten(), lat.flatten(),
                                         tlon_vec, tlat_vec, tz)
    TZ = TZv.reshape((NR, NC))

    # use good values in the new TZ in TZ0
    TZ0[~np.isnan(TZ)] = TZ[~np.isnan(TZ)]


# create a boolean mask array (True where masked = land)
# note that this is the opposite of the ROMS convention
# where mask_rho = 1. over water, and 0. over land
TZM = TZ0 > 0

if True:
    # we would also like to unmask it in the places where the
    # coastline crosses a tile, to facilitate wetting-drying
    cx = cmat['lon']
    cy = cmat['lat']
    cmask = np.isnan(cx)
    cx = cx[~cmask]
    cy = cy[~cmask]
    ii0, ii1, ifr = zfun.get_interpolant(cx, plon_vec, extrap_nan=True)
    jj0, jj1, jfr = zfun.get_interpolant(cy, plat_vec, extrap_nan=True)
    # Don't unmask extrapolated points.
    ii0 = ii0[~np.isnan(ifr) & ~np.isnan(jfr)]
    jj0 = jj0[~np.isnan(ifr) & ~np.isnan(jfr)]
    TZM[jj0, ii0] = False

#%% save the output

# save the output to NetCDF
import netCDF4 as nc
M, L = lon.shape
# get rid of old version
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist
foo = nc.Dataset(out_fn, 'w')
foo.createDimension('eta_rho', M)
foo.createDimension('xi_rho', L)
foo.createDimension('y1', M + 1)
foo.createDimension('x1', L + 1)
lat_var = foo.createVariable('lat_rho', float, ('eta_rho', 'xi_rho'))
lon_var = foo.createVariable('lon_rho', float, ('eta_rho', 'xi_rho'))
# we call these _psi_ex becasue they extend one point beyond
# that of the standard ROMS psi grid
plat_var = foo.createVariable('lat_psi_ex', float, ('y1', 'x1'))
plon_var = foo.createVariable('lon_psi_ex', float, ('y1', 'x1'))
z_var = foo.createVariable('z', float, ('eta_rho', 'xi_rho'))
mask_var = foo.createVariable('mask_rho', int, ('eta_rho', 'xi_rho'))
lat_var[:] = lat
lon_var[:] = lon
plat_var[:] = plat
plon_var[:] = plon
z_var[:] = TZ0
mask_rho = np.ones((M, L))
mask_rho[TZM == True] = 0.
mask_var[:] = mask_rho
foo.close()

#%% Test: retrieve the output

ds = nc.Dataset(out_fn)
a_z = ds.variables['z'][:]
a_mask_rho = ds.variables['mask_rho'][:]
a_lon_psi_ex = ds.variables['lon_psi_ex'][:]
a_lat_psi_ex = ds.variables['lat_psi_ex'][:]
ds.close()

# make a masked array for plotting
a_zm = ma.masked_where(a_mask_rho==0, a_z)

#%% plotting

import matplotlib.pyplot as plt

plt.close()

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
cmap1 = plt.get_cmap(name='terrain')
cs = ax.pcolormesh(a_lon_psi_ex, a_lat_psi_ex, a_zm,
                   vmin=-200, vmax=200, cmap = cmap1)
ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
zfun.dar(ax)
fig.colorbar(cs, ax=ax, extend='both')
ax.set_xlim(ax_lims[:2])
ax.set_ylim(ax_lims[-2:])
ax.set_title(gname)

plt.show()