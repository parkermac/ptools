# -*- coding: utf-8 -*-
"""
Module of functions specific to WCOFS ROMS.
"""
import netCDF4 as nc
import numpy as np
pi = np.pi


def wcofs_lonlat_2_xy(lon, lat, phi0=-57.6, theta0=37.4, regridYES = 0):

    # Regular 2D arrays of "local" coordinates x and y of WCOFS, which is 
    # a spherical grid in rotated coordinates with the pole at phi0,theta0
    # x increases along the local latitude
    # y increases in the direction opposite to local longitude
    #
    # if regridYES==1 then
    # x and y are regularized using meshgrid to eliminate roundup errors and
    # make x and y suitable to interp2. 

    phi, theta = lonlat_2_rotgrd(lon, lat, phi0, theta0)
    if regridYES:
        y, x = np.meshgrid(-phi[0,:],theta[:,0])
    else:
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

def get_basic_info(fn, only_G=False, only_S=False, only_T=False):
    """
    Gets grid, vertical coordinate, and time info from a ROMS NetCDF
    history file with full name 'fn'
    Input: the filename (with path if needed)
    Output: dicts G, S, and T
    Example calls:
    G, S, T = zfun.get_basic_info(fn)
    T = zfun.get_basic_info(fn, only_T=True)
    """
    ds = nc.Dataset(fn,'r')
    def make_G(ds):
        # get grid and bathymetry info
        g_varlist = ['h', 'x_rho', 'y_rho', 'x_u', 'y_u', 'x_v', 'y_v',
        'x_psi', 'y_psi', 'mask_rho', 'mask_u', 'mask_v', 'pm', 'pn',]
        G = dict()
        for vv in g_varlist:
            G[vv] = ds.variables[vv][:]
        G['DX'] = 1/G['pm']
        G['DY'] = 1/G['pn']
        G['M'], G['L'] = np.shape(G['x_rho']) # M = rows, L = columns
        # make the masks boolean (True = water, False = land, opposite of masked arrays!)
        G['mask_rho'] = G['mask_rho'] == 1
        G['mask_u'] = G['mask_u'] == 1
        G['mask_v'] = G['mask_v'] == 1
        return G
    def make_S(ds):
        # get vertical sigma-coordinate info (vectors are bottom to top)
        s_varlist = ['s_rho', 's_w', 'hc', 'Cs_r', 'Cs_w', 'Vtransform']
        S = dict()
        for vv in s_varlist:
            S[vv] = ds.variables[vv][:]
        S['N'] = len(S['s_rho']) # number of vertical levels
        
        # naming hack 2019.03.18
        S['sc_r'] = S['s_rho']
        S['sc_w'] = S['s_w']
        
        return S
    def make_T(ds):
        # get time info
        t_varlist = ['ocean_time', 'dstart']
        T = dict()
        for vv in t_varlist:
            T[vv] = ds.variables[vv][:]
        # find  time reference
        dstart = ds.variables['dstart']
        tu = dstart.units
        import re
        isdash = [m.start() for m in re.finditer('-', tu)]
        iscolon = [m.start() for m in re.finditer(':', tu)]
        year = int(tu[isdash[0]-4:isdash[0]])
        month = int(tu[isdash[1]-2:isdash[1]])
        day = int(tu[isdash[1]+1:isdash[1]+3])
        hour = int(tu[iscolon[0]-2:iscolon[0]])
        minute = int(tu[iscolon[1]-2:iscolon[1]])
        second = int(tu[iscolon[1]+1:iscolon[1]+3])
        import datetime
        tt = datetime.datetime(year, month, day, hour, minute, second)
        delta = datetime.timedelta(0, int(T['ocean_time']))
        T['tm0'] = tt
        T['tm'] = tt + delta
        return T
    # return results
    if only_G:
        return make_G(ds)
    elif only_S:
        return make_S(ds)
    elif only_T:
        return make_T(ds)
    else:
        return make_G(ds), make_S(ds), make_T(ds)
        