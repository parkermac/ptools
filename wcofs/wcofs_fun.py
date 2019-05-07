# -*- coding: utf-8 -*-
"""
Module of functions specific to WCOFS ROMS.
"""
import netCDF4 as nc
import numpy as np
pi = np.pi

import zfun


def wcofs_lonlat_2_xy(lon, lat, phi0=-57.6, theta0=37.4):

    # Regular 2D arrays of "local" coordinates x and y of WCOFS, which is 
    # a spherical grid in rotated coordinates with the pole at phi0,theta0
    # x increases along the local latitude
    # y increases in the direction opposite to local longitude
    #
    # 2019.05.07 TESTED against original matlab version and gives
    # identical answers

    phi, theta = lonlat_2_rotgrd(lon, lat, phi0, theta0)
    y=-phi
    x=theta
    
    return x, y

def lonlat_2_rotgrd(phi,theta,phi0,theta0):

    phi=phi*pi/180
    theta=theta*pi/180
    phi0=phi0*pi/180
    theta0=theta0*pi/180

    phi1=phi-phi0
    phi2=np.arctan2(np.sin(phi1)*np.cos(theta),
               np.cos(phi1)*np.cos(theta)*np.sin(theta0)-np.sin(theta)*np.cos(theta0))
    theta2=np.arcsin(np.cos(phi1)*np.cos(theta)*np.cos(theta0)+np.sin(theta)*np.sin(theta0))

    phi2=phi2*180/pi
    theta2=theta2*180/pi

    return phi2, theta2

def get_grid_info(fn):
    ds = nc.Dataset(fn,'r')
    # get grid and bathymetry info
    g_varlist = ['h', 'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v',
    'lon_psi', 'lat_psi', 'mask_rho', 'mask_u', 'mask_v', 'mask_psi',
    'pm', 'pn', 'angle']
    G = dict()
    for vv in g_varlist:
        G[vv] = ds.variables[vv][:]
    G['DX'] = 1/G['pm']
    G['DY'] = 1/G['pn']
    G['M'], G['L'] = np.shape(G['lon_rho']) # M = rows, L = columns
    # make the masks boolean (True = water, False = land, opposite of masked arrays!)
    G['mask_rho'] = G['mask_rho'] == 1
    G['mask_u'] = G['mask_u'] == 1
    G['mask_v'] = G['mask_v'] == 1
    ds.close()
    
    return G
    
def get_itp_dict(sta_dict, sta_list, GG, Gcut):
    # get interpolants
    itp_dict = dict()
    for sta_name in sta_list:#sta_dict.keys():
        # target position
        Lon = sta_dict[sta_name][0]
        Lat = sta_dict[sta_name][1]
        [X, Y]=wcofs_lonlat_2_xy(Lon,Lat);
        X = np.array([X])
        Y = np.array([Y])
        # get interpolants for this point
        Xi0 = dict(); Yi0 = dict()
        Xi1 = dict(); Yi1 = dict()
        Aix = dict(); Aiy = dict()
        for grd in ['rho', 'u', 'v']:
            xx = GG['x_' + grd][0,:]
            yy = GG['y_' + grd][:,0]
            mm = Gcut['mask_' + grd]
            xi0, xi1, xfr = zfun.get_interpolant(X, xx, extrap_nan=True)
            yi0, yi1, yfr = zfun.get_interpolant(Y, yy, extrap_nan=True)
        
            # check for masked points
            # Note: becasue we do this separately for each grid we can end up with
            # a situation where the rho, u, and v points are not at the same location.
            # for now we will just accept this because it only applies to a few shallow moorings.
            for nudges in range(5):
                mask_test = mm[yi0[0]:yi1[0]+1, xi0[0]:xi1[0]+1].flatten()
                if (mask_test == False).any():
                    print(' Warning - masked points for ' + sta_name + ': grid = ' + grd)
                    print(' nudge to the west # ' + str(nudges+1))
                    xi0 -= 1
                    xi1 -= 1
            
            Xi0[grd] = xi0
            Yi0[grd] = yi0
            Xi1[grd] = xi1
            Yi1[grd] = yi1
            # create little arrays that are used in the actual interpolation
            Aix[grd] = np.array([1-xfr, xfr]).reshape((1,1,2))
            Aiy[grd] = np.array([1-yfr, yfr]).reshape((1,2))
        itp_dict[sta_name] = (Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
    return itp_dict
    
def get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy):
    # convenience function used during the extraction
    dims = ds[vv].dimensions
    if 'eta_rho' in dims:
        grd = 'rho'
    elif 'eta_u' in dims:
        grd = 'u'
    elif 'eta_v' in dims:
        grd = 'v'
    else:
        print('grid error!')
    xi0 = Xi0[grd]; yi0 = Yi0[grd]
    xi1 = Xi1[grd]; yi1 = Yi1[grd]
    aix = Aix[grd]; aiy = Aiy[grd]
    xi01 = np.array([xi0, xi1]).flatten()
    yi01 = np.array([yi0, yi1]).flatten()
    return xi01, yi01, aix, aiy
