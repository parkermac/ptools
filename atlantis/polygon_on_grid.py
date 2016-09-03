# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 13:02:17 2016

@author: PM5

Code to create indices of grid cells that best approximate a polygon.

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()
import zfun
import zrfun

pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

import pandas as pd
import matplotlib.pyplot as plt
#import cmocean as cmo
import netCDF4 as nc4
import numpy as np

from datetime import datetime

#%% get history file

in_dir0 = Ldir['parent'] + 'roms/output/salish_2006_4_lp/'
dt = datetime(2006,7,29)
f_string = 'f' + dt.strftime('%Y.%m.%d')
in_dir = in_dir0 + f_string + '/'
fn = in_dir + 'low_passed.nc'

[G] = zrfun.get_basic_info(fn, getS=False, getT=False)
ds = nc4.Dataset(fn)

lon = G['lon_psi']
lat = G['lat_psi']
Lon = lon[0,:].flatten()
Lat = lat[:,0].flatten()

#%% load Atlantis polygon info

pfn = (Ldir['parent'] + 'PROJECTS/LLTK/Atlantis/Puget_Sound_HydroAtlantis/' +
        'AtlantisBoxInfo_toParker.xlsx')

df = pd.read_excel(pfn, sheetname='BoxVertices')

#%% get a polygon

npoly = 0 # go to df.box_id.max()

a = df[df.box_id==npoly]
lon_poly = a.Long.values
lat_poly = a.Lat.values

# initialize result dicts
glon_dict = dict()
glat_dict = dict()

def ll2xy(lon, lat, lon0, lat0, R, clat):
    # should work for lon, lat scalars or arrays
    x = R * clat * np.pi * (lon - lon0) / 180
    y = R * np.pi * (lat - lat0) / 180
    return x, y

def get_dist_normal(x1, x2, y1, y2, xp1, yp1):
    # find shortest distance from xp1, yp1 to the segment

    if x2 < x1:
        (x1, x2) = (x2, x1)
        (y1, y2) = (y2, y1)

    dx = x2-x1
    dy = y2-y1
    if dy == 0:
        xp2 = xp1
        yp2 = y1
    elif dx == 0:
        yp2 = yp1
        xp2 = x1
    else:
        m = dy/dx
        b = y1 - m*x1
        bp = yp1 + (1/m)*xp1
        xp2 = (bp - b)/(m + (1/m))
        yp2 = m*xp2 + b
    dxp = xp2 - xp1
    dyp = yp2 - yp1
    dist = np.sqrt(dxp**2 + dyp**2)
    return dist

EW = np.array([0, 1, 0, -1], dtype=int)
NS = np.array([1, 0, -1, 0], dtype=int)

for iseg in range(len(lon_poly) - 1):
    print('**************** segment %d ************************' % (iseg))


    lon0 = lon_poly[iseg]
    lat0 = lat_poly[iseg]

    # convert to meters from the first point in the polygon
    R = zfun.earth_rad(lat0)
    clat = np.cos(np.pi*lat0/180)

    # the x, y positions of the polygon vertices
    x, y = ll2xy(lon_poly, lat_poly, lon0, lat0, R, clat)

    # get index of nearest point in the psi grid
    ix0 = zfun.find_nearest_ind(Lon, lon0)
    iy0 = zfun.find_nearest_ind(Lat, lat0)

    # x, y position of the starting point on the psi grid
    xp1, yp1 = ll2xy(Lon[ix0], Lat[iy0], lon0, lat0, R, clat)

    # initialize result vectors
    glon_dict[iseg] = [Lon[ix0]]
    glat_dict[iseg] = [Lat[iy0]]

    x1 = x[iseg]
    x2 = x[iseg+1]
    y1 = y[iseg]
    y2 = y[iseg+1]

    Dist = np.nan * np.ones(4)
    DistN = np.nan * np.ones(4)

    do_next = False
    while do_next == False:

        xp0, yp0 = ll2xy(Lon[ix0], Lat[iy0], lon0, lat0, R, clat)
        Dist0 = np.sqrt((x2-xp0)**2 + (y2-yp0)**2)

        # find the distance of the four neighboring points to the segment
        # arranged N, E, S, W
        for jj in range(4):
            xp1, yp1 = ll2xy(Lon[ix0+EW[jj]], Lat[iy0+NS[jj]], lon0, lat0, R, clat)
            Dist[jj] = np.sqrt((x2-xp1)**2 + (y2-yp1)**2)
            DistN[jj] = get_dist_normal(x1, x2, y1, y2, xp1, yp1)
        dDist = Dist - Dist0

        mask = dDist < 0 # mask True for closer points

        if sum(mask) == 0:
            do_next = True
            break

        DistN[~mask] = np.nan
        ii_nesw = np.nanargmin(DistN)

        ix_next = ix0+EW[ii_nesw]
        iy_next = iy0+NS[ii_nesw]

        xp_next, yp_next = ll2xy(Lon[ix_next], Lat[iy_next], lon0, lat0, R, clat)
        Dist_next = np.sqrt((x2-xp_next)**2 + (y2-yp_next)**2)

        print('Dist0 = ' + str(Dist0))
        print('Dist_next = ' + str(Dist_next))
        print('Dist: ' + str(Dist))
        print('DistN: ' + str(DistN))
        print(' ')

        if (Dist_next < Dist0):
            ix0 = ix_next
            iy0 = iy_next
            glon_dict[iseg].append(Lon[ix0])
            glat_dict[iseg].append(Lat[iy0])
        else:
            do_next = True

#%% plotting

if True:

    plt.close('all')

    fig = plt.figure(figsize=(12, 14))
    ax = fig.add_subplot(111)
    #ax.pcolormesh(lon, lat, G['h'][1:-1, 1:-1], cmap=cmo.cm.deep, alpha=.2)
    pfun.add_coast(ax)
    ax.axis([-124, -122, 46.8, 49.2])
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

    ax.plot(lon_poly, lat_poly,'-*r')
    ax.plot(lon_poly[0], lat_poly[0],'-or', markersize=12)
    ax.plot(lon_poly[:2], lat_poly[:2],'-g', linewidth=2)

    for iseg in glon_dict.keys():
        ax.plot(glon_dict[iseg], glat_dict[iseg], '-ob')
        ax.plot(glon_dict[iseg][0], glat_dict[iseg][0], '-ob', markersize=12, alpha=.2)
        ax.plot(glon_dict[iseg][-1], glat_dict[iseg][-1], '-*y', markersize=12)

    V = np.ones((len(lon_poly),2))
    V[:,0] = lon_poly
    V[:,1] = lat_poly
    import matplotlib
    P = matplotlib.path.Path(V)
    Rlon = G['lon_rho'].flatten()
    Rlat = G['lat_rho'].flatten()
    R = np.ones((len(Rlon),2))
    R[:,0] = Rlon
    R[:,1] = Rlat
    RR = P.contains_points(R)
    ax.plot(Rlon[RR], Rlat[RR], 'gs')

    plt.show()