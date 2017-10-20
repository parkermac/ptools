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
import gfun
reload(gfun)
Gr =gfun.gstart()
import gfun_utility as gfu
reload(gfu)

import numpy as np
import pickle

import Lfun
import zfun

#%%
Lfun.make_dir(Gr['gdir'], clean=True)

fn = 'grid_m00_r00_s00_x00.nc'
out_fn = Gr['gdir'] + fn
print(50*'*')
print(out_fn)

dch = gfun.default_choices(Gr)

#%% GRID DEFINITIONS

# vectors to define the plaid grid
# start with cell corners (like an extended psi grid)

if Gr['gridname'] == 'big1':
    # the next evolution: big domain, high-res
    dch =  gfun.default_choices(Gr, wet_dry=True) # override for wet_dry
    dch['z_land'] = 2
    maxres = 4000
    medres = 1000
    minres = 300
    lon_list = [-130, -124.5, -123.5, -122]
    x_res_list = [maxres, medres, minres, minres]
    lat_list = [42, 46, 47, 49, 50, 52]
    y_res_list = [maxres, medres, minres, minres, medres, maxres]
    plon_vec, plat_vec = gfu.stretched_grid(lon_list, x_res_list,
                                        lat_list, y_res_list)
    dch['nudging_edges'] = ['north', 'south', 'west']

if Gr['gridname'] == 'big2':
    # the next-next evolution: big domain, med-high-res
    dch =  gfun.default_choices(Gr, wet_dry=True) # override for wet_dry
    dch['z_land'] = 2
    dch['min_depth'] = -5 # negative means land
    maxres = 4000
    medres = 1000
    minres = 600
    lon_list = [-130, -124.5, -123.5, -122]
    x_res_list = [maxres, medres, minres, minres]
    lat_list = [42, 46, 47, 49, 50, 52]
    y_res_list = [maxres, medres, minres, minres, medres, maxres]
    plon_vec, plat_vec = gfu.stretched_grid(lon_list, x_res_list,
                                        lat_list, y_res_list)
    dch['nudging_edges'] = ['north', 'south', 'west']

    
elif Gr['gridname'] == 'sal0':
    # start of a salish nest grid
    aa = [-124, -122, 47, 49]
    res = 300 # target resolution (m)
    plon_vec, plat_vec = gfu.simple_grid(aa, res)
    dch['t_list'] = ['cascadia/cascadia_gridded.nc',
             'psdem/PS_183m.nc',
             'ttp_patch/TTP_Regional_27m_patch.nc']
    dch['nudging_edges'] = ['north', 'west']
    
elif Gr['gridname'] == 'hc0':
    # mid Hood Canal nest
    aa = [-123, -122.55, 47.5, 47.9]
    res = 50 # target resolution (m)
    plon_vec, plat_vec = gfu.simple_grid(aa, res)
    dch['t_list'] = ['psdem/PS_27m.nc']
    dch['nudging_edges'] = ['north', 'west']
             
elif Gr['gridname'] == 'cas1': # An extended version of cascadia1
    maxres = 5000
    medres = 3000
    minres = 1500
    lon_list = [-127.4, -126, -124, -122]
    x_res_list = [maxres, medres, minres, minres]
    lat_list = [42, 47.9, 49, 50]
    y_res_list = [medres, minres, minres, medres]
    plon_vec, plat_vec = gfu.stretched_grid(lon_list, x_res_list,
                                        lat_list, y_res_list)
    dch['use_cell_average'] = True
    dch['nudging_edges'] = ['south', 'west']
                                        
elif Gr['gridname'] == 'aestus1': # idealized model
    lon_list = [-1, 0, 1, 2, 3]
    x_res_list = [5000, 1000, 1000, 5000, 5000]
    lat_list = [44, 44.9, 45.1, 46]
    y_res_list = [5000, 1000, 1000, 5000]
    plon_vec, plat_vec = gfu.stretched_grid(lon_list, x_res_list,
                                        lat_list, y_res_list)
    dch['analytical'] = True
    dch['nudging_edges'] = ['north', 'south', 'west']

# save the default choices for use by other code
pickle.dump(dch, open(Gr['gdir'] + 'choices.p', 'wb'))

plon, plat = np.meshgrid(plon_vec, plat_vec)
ax_lims = (plon_vec[0], plon_vec[-1], plat_vec[0], plat_vec[-1])

# make box centers
lon_vec = plon_vec[:-1] + np.diff(plon_vec)/2
lat_vec = plat_vec[:-1] + np.diff(plat_vec)/2
lon, lat = np.meshgrid(lon_vec, lat_vec)
NR, NC = lon.shape

# initialize the final bathymetry array
z = np.nan * lon
if dch['analytical']==True:
    if Gr['gridname'] == 'aestus1':
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
    if dch['do_cell_average']:
        # m is the start of a mask: 1=water, 0=land
        m = np.nan * lon
        for t_file in dch['t_list']:
            t_fn = dch['t_dir'] + t_file
            print('\nOPENING BATHY FILE: ' + t_file)
            tlon_vec, tlat_vec, tz = gfu.load_bathy_nc(t_fn)
            if isinstance(tz, np.ma.masked_array):
                tz1 = tz.data
                tz1[tz.mask==True] = np.nan
                tz = tz1
            # Apply the offset here instead of at the end because it matters
            # for the estimated mask field.
            if dch['use_z_offset']:
                tz = tz + dch['z_offset']
            tm = np.ones_like(tz)
            tm[tz>0] = 0.
            # average in grid cells
            xi0, xi1, xf = zfun.get_interpolant(tlon_vec,plon_vec, extrap_nan=True)
            yi0, yi1, yf = zfun.get_interpolant(tlat_vec,plat_vec, extrap_nan=True)
            z_part = np.nan * np.ones((NR,NC))
            m_part = np.nan * np.ones((NR,NC))
            tNR, tNC = tz.shape
            itx = np.arange(tNC)
            jty = np.arange(tNR)
            for ii in range(NC):
                for jj in range(NR):
                    ix = itx[xi0==ii]
                    jy = jty[yi0==jj]
                    if ix.size>0 and jy.size>0:
                        z_part[jj, ii] = np.nanmean(tz[jy[0]:jy[-1], ix[0]:ix[-1]])
                        m_part[jj, ii] = np.nanmean(tm[jy[0]:jy[-1], ix[0]:ix[-1]])
                    else:
                        pass
            # put good values of z_part in z
            z[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]
            m[~np.isnan(m_part)] = m_part[~np.isnan(m_part)]            
    else:
        # m is the start of a mask: 1=water, 0=land
        m = np.ones_like(lon)
        for t_file in dch['t_list']:
            t_fn = dch['t_dir'] + t_file
            print('\nOPENING BATHY FILE: ' + t_file)
            tlon_vec, tlat_vec, tz = gfu.load_bathy_nc(t_fn)
            if isinstance(tz, np.ma.masked_array):
                tz1 = tz.data
                tz1[tz.mask==True] = np.nan
                tz = tz1
            tlon, tlat = np.meshgrid(tlon_vec, tlat_vec)
            z_part = zfun.interp2(lon, lat, tlon, tlat, tz)
            # put good values of z_part in z
            z[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]
        if dch['use_z_offset']:
            z = z + dch['z_offset']

#%% save the output to NetCDF
gfu.make_nc(out_fn, plon, plat, lon, lat, z, m, dch)



