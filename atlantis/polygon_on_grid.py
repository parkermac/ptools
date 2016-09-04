# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 13:02:17 2016

@author: PM5

Code to create indices of grid cells that best approximate a polygon. Created
for the Atlantis project, funded by LLTK.

RESULT: this does what appears to be a perfect job of defining the psi-grid
points of each segment, which share beginning and ending points, and the
rho-grid points inside the original polygon are also exactly those
contained by the psi-grid segments.  Nice!

"""

#%% Imports

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
import matplotlib.path as mpath
import numpy as np
from datetime import datetime

#%% get ROMS history file

in_dir0 = Ldir['parent'] + 'roms/output/salish_2006_4_lp/'
dt = datetime(2006,7,29)
f_string = 'f' + dt.strftime('%Y.%m.%d')
in_dir = in_dir0 + f_string + '/'
fn = in_dir + 'low_passed.nc'
[G] = zrfun.get_basic_info(fn, getS=False, getT=False)
lon = G['lon_psi']
lat = G['lat_psi']
# get the vectors that describe the model plaid grid
Lon = lon[0,:].flatten()
Lat = lat[:,0].flatten()

#%% load Atlantis polygon info into a DataFrame

pfn = (Ldir['parent'] + 'PROJECTS/LLTK/Atlantis/Puget_Sound_HydroAtlantis/' +
        'AtlantisBoxInfo_toParker.xlsx')
df = pd.read_excel(pfn, sheetname='BoxVertices')

#%% get a polygon

# to iterate over all polygons in the DataFrame, go to df.box_id.max()
npoly = 0
this_poly = df[df.box_id==npoly]
lon_poly = this_poly.Long.values
lat_poly = this_poly.Lat.values

# Initialize result dicts
#
# Each item is a list of lat or lon values on the ROMS psi grid for a
# single segment of the polygon.
# (eventually we want to save indices as well)
glon_dict = dict()
glat_dict = dict()

def ll2xy(lon, lat, lon0, lat0, R, clat):
    # This converts lon, lat into meters relative to lon0, lat0.
    # It should work for lon, lat scalars or arrays.
    x = R * clat * np.pi * (lon - lon0) / 180
    y = R * np.pi * (lat - lat0) / 180
    return x, y

def get_dist_normal(x1, x2, y1, y2, xp1, yp1):
    # This finds the shortest distance from xp1, yp1 to the segment
    # defined by points (x1, y1) and (x2, y2).
    if x2 < x1:
        # make sure x is increasing (inside the function)
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

# arrays to define which way to increase indices, arranged in order
# NESW
EW = np.array([0, 1, 0, -1], dtype=int)
NS = np.array([1, 0, -1, 0], dtype=int)

for iseg in range(len(lon_poly) - 1):
    # iterate over all segments in the polygon
    #print('**************** segment %d ************************' % (iseg))
    lon0 = lon_poly[iseg]
    lat0 = lat_poly[iseg]
    # convert to meters from the first point in the polygon
    R = zfun.earth_rad(lat0)
    clat = np.cos(np.pi*lat0/180)
    # the x, y positions of all the polygon vertices
    # (although we only use two for each segment)
    x, y = ll2xy(lon_poly, lat_poly, lon0, lat0, R, clat)
    # get index of nearest point on the psi grid to lon0, lat0
    ix0 = zfun.find_nearest_ind(Lon, lon0)
    iy0 = zfun.find_nearest_ind(Lat, lat0)
    # x, y position of the starting point on the psi grid
    xp1, yp1 = ll2xy(Lon[ix0], Lat[iy0], lon0, lat0, R, clat)

    # load first point in result lists
    glon_dict[iseg] = [Lon[ix0]]
    glat_dict[iseg] = [Lat[iy0]]
    # get the endpoints of this segment
    x1 = x[iseg]
    x2 = x[iseg+1]
    y1 = y[iseg]
    y2 = y[iseg+1]
    # Initialize the array of the distances of the four points surrounding
    # the current psi grid point relative to the end point of the segment.
    # This gives direction to the path.
    Dist = np.nan * np.ones(4)
    # Initialize the array of the distances of the four points surrounding
    # the current psi grid point relative to the segment.
    # This keeps the path close to the segment.
    DistN = np.nan * np.ones(4)

    # iterate on this segment until we get to where the distance to the
    # segment endpoint is no longer decreasing
    do_next_segment = False
    while do_next_segment == False:
        # find the distance of the current point to to segment endpoint
        xp0, yp0 = ll2xy(Lon[ix0], Lat[iy0], lon0, lat0, R, clat)
        Dist0 = np.sqrt((x2-xp0)**2 + (y2-yp0)**2)
        # Find the distance of the four neighboring points to the endpoint
        # and to the segment, arranged in order N, E, S, W.
        for jj in range(4):
            xp1, yp1 = ll2xy(Lon[ix0+EW[jj]], Lat[iy0+NS[jj]], lon0, lat0, R, clat)
            Dist[jj] = np.sqrt((x2-xp1)**2 + (y2-yp1)**2)
            DistN[jj] = get_dist_normal(x1, x2, y1, y2, xp1, yp1)
        # only work with distance in the 4 choices which are closer than
        # the current center point
        dDist = Dist - Dist0
        mask = dDist < 0 # mask is True for closer points
        # If no points are closer then we are done with this segment
        if sum(mask) == 0:
            do_next = True
            break
        # otherwise we first mask out points which are not closer to the endpoint
        DistN[~mask] = np.nan
        # and then choose the one that is closest to the segment
        ii_nesw = np.nanargmin(DistN)
        # and use this to figure out the index of the next point on the path
        ix_next = ix0+EW[ii_nesw]
        iy_next = iy0+NS[ii_nesw]
        xp_next, yp_next = ll2xy(Lon[ix_next], Lat[iy_next], lon0, lat0, R, clat)
        # this is the distance from the next point to the end
        Dist_next = np.sqrt((x2-xp_next)**2 + (y2-yp_next)**2)
        if (Dist_next < Dist0):
            # If this point has made progress toward the goal,
            # then append it to the results lists, and update ix0 and iy0.
            ix0 = ix_next
            iy0 = iy_next
            glon_dict[iseg].append(Lon[ix0])
            glat_dict[iseg].append(Lat[iy0])
        else:
            # otherwise, we are done with this segment
            do_next_segment = True

# Trim extra shared points at the end or biginning of each segment
for iseg in glon_dict.keys():
    xx0 = glon_dict[iseg]
    yy0 = glat_dict[iseg]
    if iseg == len(glon_dict.keys()) - 1:
        iseg1 = 0
    else:
        iseg1 = iseg+1
    xx1 = glon_dict[iseg1]
    yy1 = glat_dict[iseg1]
    keep_trimming = True
    while keep_trimming == True:
        A0 = (xx0[-2], xx0[-1], yy0[-2], yy0[-1])
        A1 = (xx1[1], xx1[0], yy1[1], yy1[0])
        if A0 == A1:
            xx0.pop()
            yy0.pop()
            xx1.pop(0)
            yy1.pop(0)
        else:
            keep_trimming = False
    glon_dict[iseg] = xx0
    glat_dict[iseg] = yy0
    glon_dict[iseg1] = xx1
    glat_dict[iseg1] = yy1

# find points on the rho grid that are inside the polygon
V = np.ones((len(lon_poly),2))
V[:,0] = lon_poly
V[:,1] = lat_poly
P = mpath.Path(V)
Rlon = G['lon_rho'].flatten()
Rlat = G['lat_rho'].flatten()
R = np.ones((len(Rlon),2))
R[:,0] = Rlon
R[:,1] = Rlat
RR = P.contains_points(R) # boolean

#%% plotting

if True:

    plt.close('all')

    fig = plt.figure(figsize=(12, 14))
    ax = fig.add_subplot(111)
    pfun.add_coast(ax)
    ax.axis([-124, -122, 46.8, 49.2])
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

    # plot the original polygon, with a star at the start and
    # a segment indicating direction
    ax.plot(lon_poly, lat_poly,'-r')
    ax.plot(lon_poly[0], lat_poly[0],'-*m', markersize=12)
    ax.plot(lon_poly[:2], lat_poly[:2],'-m', linewidth=3)

    for iseg in glon_dict.keys():
        # plot the points on the psi grid for each segment (stairstep lines)
        ax.plot(glon_dict[iseg], glat_dict[iseg], '-xb')
        # starting points
        ax.plot(glon_dict[iseg][0], glat_dict[iseg][0], '-*b', markersize=12)
        # ending points
        ax.plot(glon_dict[iseg][-1], glat_dict[iseg][-1], '-oy', markersize=12, alpha=.2)

    # plot points on the rho grid that are inside the polygon
    ax.plot(Rlon[RR], Rlat[RR], 'go')

    plt.show()