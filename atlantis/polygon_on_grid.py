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

It also makes a dict Shared_faces" that tells which face of which polygon
each face of a given polygon is connected to.

Reworked 11/12/2016 to omit clipping, and try a different way to determine
the direction of "in" so that all polygons with shared faces have exactly
matching info.

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

import pandas as pd
import matplotlib.path as mpath
import numpy as np
from datetime import datetime
import pickle

#%% setup output location

whichyear = 2005
if Ldir['env'] == 'pm_mac': # mac version
    if whichyear == 2006:
        R_in_dir0 = Ldir['parent'] + 'roms/output/salish_2006_4_lp/'
        out_dir0 = Ldir['parent'] + 'ptools_output/atlantis_mac_2006/'
    elif whichyear == 2005:
        R_in_dir0 = Ldir['parent'] + 'roms/output/salish_2005_1_lp/'
        out_dir0 = Ldir['parent'] + 'ptools_output/atlantis_mac_2005/'
elif Ldir['env'] == 'pm_fjord': # fjord version
    if whichyear == 2006:
        R_in_dir0 = '/boildat1/parker/roms/output/salish_2006_4_lp/'
        out_dir0 = Ldir['parent'] + 'ptools_output/atlantis_fjord_2006/'
    elif whichyear == 2005:
        R_in_dir0 = '/boildat1/parker/roms/output/salish_2005_1_lp/'
        out_dir0 = Ldir['parent'] + 'ptools_output/atlantis_fjord_2005/'
Lfun.make_dir(out_dir0)
out_dir = out_dir0 + 'gridded_polygons/'
Lfun.make_dir(out_dir, clean=True)

#%% get ROMS history file

dt = datetime(whichyear,7,29)
f_string = 'f' + dt.strftime('%Y.%m.%d')
R_in_dir = R_in_dir0 + f_string + '/'
R_fn = R_in_dir + 'low_passed.nc'
G = zrfun.get_basic_info(R_fn, only_G=True)
lon = G['lon_psi']
lat = G['lat_psi']
# get the vectors that describe the model plaid grid
Lon = lon[0,:].flatten()
Lat = lat[:,0].flatten()

#%% load Atlantis polygon info into a DataFrame

pfn = (Ldir['parent'] + 'ptools_data/atlantis/' +
        'AtlantisBoxInfo_toParker.xlsx')
df = pd.read_excel(pfn, sheetname='BoxVertices')
NPOLY = df.box_id.max() + 1

#%% define some functions

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

#%% process one or more polygons

gpoly_dict = dict()

for npoly in range(NPOLY):

    print('npoly = ' + str(npoly))

    this_poly = df[df.box_id==npoly]
    lon_poly = this_poly.Long.values
    lat_poly = this_poly.Lat.values

    # Initialize result dicts
    #
    # Each item is a list indices into Lon or Lat on the ROMS psi grid for a
    # single segment of the polygon.
    # i = column, j = row
    i_dict = dict()
    j_dict = dict()

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

        # save first point in result lists
        i_dict[iseg] = [ix0]
        j_dict[iseg] = [iy0]
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
                i_dict[iseg].append(ix0)
                j_dict[iseg].append(iy0)
            else:
                # otherwise, we are done with this segment
                do_next_segment = True

    # find points on the rho grid that are inside the new polygon
    #
    # if we constructed V this way is gives almost the same results
    # but there can be a few inconsistencies
    # V = np.ones((len(lon_poly),2))
    # V[:,0] = lon_poly
    # V[:,1] = lat_poly
    #
    # doing it with the stair step psi-grid points means there are no
    # inconsistencies
    plon_poly = np.array([])
    plat_poly = np.array([])
    for iseg in i_dict.keys():
        plon_poly = np.concatenate((plon_poly, Lon[i_dict[iseg]]))
        plat_poly = np.concatenate((plat_poly, Lat[j_dict[iseg]]))
    V = np.ones((len(plon_poly),2))
    V[:,0] = plon_poly
    V[:,1] = plat_poly
    P = mpath.Path(V)
    Rlon = G['lon_rho'].flatten()
    Rlat = G['lat_rho'].flatten()
    R = np.ones((len(Rlon),2))
    R[:,0] = Rlon
    R[:,1] = Rlat
    RR = P.contains_points(R) # boolean
    # create arrays of i (column) and j (row) indices
    i_rho = np.arange(G['L']).reshape((1,G['L'])).repeat(G['M'], axis=0)
    j_rho = np.arange(G['M']).reshape((G['M'],1)).repeat(G['L'], axis=1)
    # pack indices that are inside the polygon
    # as a numpy int array, with two columns, packed in order j,i
    ji_rho_in = np.array([j_rho.flatten()[RR], i_rho.flatten()[RR]],
                         dtype=int).T

    # pack perimiter data into a dict of arrays, one array for each segment
    # and each row in an (int) array is:
    # [j index in grid, i index in grid, 0=u-grid or 1=v-grid, +/-1 sign]
    # where the sign is positive if positive velocity at that point
    # would be into the polygon (as determined by ji_rho_in).
    per_dict = dict()
    for iseg in i_dict.keys():
        iis = i_dict[iseg]
        jjs = j_dict[iseg]
        
        # determine the slope of the line defining this face
        dx = iis[-1] - iis[0]
        dy = jjs[-1] - jjs[0]
        # and assume we are moving CCW around each polygon so that
        # "in" is to the left.
        
        per_dict[iseg] = np.zeros((len(iis)-1,4), dtype=int)
        # if the segment is a single psi point then we end up with
        # per_dict[iseg] = array([], shape=(0, 4), dtype=int64)
        for nn in range(len(iis)-1):
            # get indices for two points on the psi-grid
            ip0 = iis[nn]
            ip1 = iis[nn+1]
            jp0 = jjs[nn]
            jp1 = jjs[nn+1]
            # determine if they encompass a u or v point and assign j and i
            if ip0 == ip1:
                uv = 0 # u-grid
                this_i = ip0
                this_j = np.max([jp0, jp1])
            elif jp0 == jp1:
                uv = 1 # v_grid
                this_i = np.max([ip0, ip1])
                this_j = jp0
            else:
                print('psi index error')
            # determine sign for flow into the polygon by looking for
            # points on the rho grid that are inside
            if uv == 0: # u-grid
                if dy > 0:
                    pm = -1
                elif dy < 0:
                    pm = 1
                else:
                    print('\nnpoly=%d nface=%d' % (npoly, iseg))
                    print('found u-grid point where dy = 0')
            elif uv == 1: # v-grid
                if dx > 0:
                    pm = 1
                elif dx < 0:
                    pm = -1
                else:
                    print('\nnpoly=%d nface=%d' % (npoly, iseg))
                    print('found v-grid point where dx = 0')
                
            per_dict[iseg][nn, :] = [this_j, this_i, uv, pm]

    # save output
    gpoly_dict[npoly] = {'per_dict':per_dict, 'ji_rho_in':ji_rho_in,
            'i_dict':i_dict, 'j_dict':j_dict}

# write output to disk
pickle.dump(gpoly_dict, open(out_dir + 'gpoly_dict.p', 'wb'))

#%% *** shared faces code ***

# First, make a dict of all the faces coordinates
poly_face = dict()
# dict structure
# (npoly, nface) : (lon1, lon0, lat1,lat0)
for npoly in range(NPOLY):
    this_poly = df[df.box_id==npoly]
    lon_poly = this_poly.Long.values
    lat_poly = this_poly.Lat.values
    for nface in range(len(lon_poly) - 1):
        # Note that here we take the segment (face) ends in reverse order,
        # so that they will match with those taken in regular order
        # in the later loop (*).
        poly_face[(npoly, nface)] = (lon_poly[nface+1], lon_poly[nface],
            lat_poly[nface+1], lat_poly[nface]) 

# Next initialize a dict to relate each polygon and face
# to the one it connects with.
shared_faces = dict()
# (npoly, nface): (npoly: nface)

# Finally find the shared faces and write them to the dict.
for npoly in range(NPOLY):
    this_poly = df[df.box_id==npoly]
    lon_poly = this_poly.Long.values
    lat_poly = this_poly.Lat.values    
    # figure out which faces of which polygons this one connects to
    this_poly_face = dict()
    for nface in range(len(lon_poly) - 1): # (*)
        this_poly_face[(npoly, nface)] = (lon_poly[nface], lon_poly[nface+1],
            lat_poly[nface], lat_poly[nface+1])        
    for this_pf in this_poly_face.keys():
        for pf in poly_face.keys():
            if (this_poly_face[this_pf] == poly_face[pf]):
                shared_faces[this_pf] = pf
                
# write output to disk
pickle.dump(shared_faces, open(out_dir + 'shared_faces.p', 'wb'))


