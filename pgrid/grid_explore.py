# -*- coding: utf-8 -*-
"""
Tool for interactively exploring the bathymetry of a
grid file.

Can be run im ipython with a user-specified grid file

run grid_explore.py -g sal0


"""

import gfun
from importlib import reload
reload(gfun)
Gr =gfun.gstart()
# running gfun.gstart() sets the path to include pfun and zfun
import gfun_utility as gfu
reload(gfu)
import gfun_plotting as gfp
reload(gfp)
import pfun
import zfun
reload(zfun)

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import os
import shutil

import Lfun
Ldir = Lfun.Lstart()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', nargs='?', default='',
        type=str)
args = parser.parse_args()
if len(args.gridname) > 0:
    Gr = gfun.gstart(gridname=args.gridname)
else:
    Gr = gfun.gstart()
    
# select grid file
using_old_grid = False
# Set this to True to look at grids we have already created,
# e.g. ones currently in use for LiveOcean.
# Set it to False when interacting with grids from pgrid_output.
if using_old_grid==True:
    fn = gfun.select_file(Gr, using_old_grid=True)
    in_fn = fn
elif using_old_grid==False:
    fn = gfun.select_file(Gr)
    in_fn = Gr['gdir'] + fn

# get fields
ds = nc.Dataset(in_fn)
H = ds.variables['h'][:]
mask_rho = ds.variables['mask_rho'][:]
plon = ds.variables['lon_psi_ex'][:]
plat = ds.variables['lat_psi_ex'][:]
lon = ds.variables['lon_rho'][:]
lat = ds.variables['lat_rho'][:]
DA = (1/ds['pm'][:]) * (1/ds['pn'][:])
DA[mask_rho==0] = np.nan
H[mask_rho==0] = np.nan
ds.close()

# get distances
XM, YM = zfun.ll2xy(lon, lat, np.mean(lon[0,:]), np.mean(lat[:,0]))

# flip to work with imshow
h = np.flipud(H)
da = np.flipud(DA)
xm = np.flipud(XM)
ym = np.flipud(YM)
m = np.flipud(mask_rho) # mask_rho: 1 = water, 0 = land

NR, NC = h.shape

# PLOTTING

# set up the axes
plt.close('all')
fig = plt.figure(figsize=(13,8)) # (13,8) is good for my laptop
ax1 = plt.subplot2grid((1,3), (0,0), colspan=2) # map
ax2 = plt.subplot2grid((1,3), (0,2), colspan=1) # buttons

#%% initialize the data plot

cmap1 = plt.get_cmap(name='rainbow_r') # terrain
tvmin = -20
tvmax = 200
cs = ax1.imshow(h, interpolation='nearest', vmin=tvmin, vmax=tvmax, cmap = cmap1)
fig.colorbar(cs, ax=ax1, extend='both')
aa = ax1.axis()

# add the coastline
clon, clat = pfun.get_coast()
cx0, cx1, cxf = zfun.get_interpolant(clon, lon[0,:], extrap_nan=True)
cy0, cy1, cyf = zfun.get_interpolant(clat, lat[:,0], extrap_nan=True)
ax1.plot(cx0 + cxf, NR - (cy0 + cyf) - 1, '-k')

# add rivers
gfp.edit_mask_river_tracks(Gr, NR, ax1)

ax1.axis(aa)

# create control buttons
# list is organized from bottom to top
blist = ['start', 'pause','polyInfo', 'lineLength', 'startPoly',
         'done']
# nicer names
Blist = ['Start', 'Pause','Polygon Info',
         'Line Length', 'Start Polygon\nLine',
         'Done']
NB = len(blist) # number of buttons
ybc = np.arange(NB+1) - .5
offset = 1e5 # kludgey way to distinguish buttons from topography
xbc = np.arange(2) + offset
XBC, YBC = np.meshgrid(xbc, ybc)
ZB = np.arange(1,NB+1).reshape(NB,1)
cmap2 = plt.get_cmap(name='viridis')
ax2.pcolormesh(XBC,YBC,ZB, vmin=0, vmax=NB, cmap = cmap2)
ax2.set_axis_off()

def addButtonLabel(ax, plon, plat, nb, lab, tcol='k'):
    # draw and label buttons
    pad = .1
    ax.add_patch(
        plt.Rectangle((plon[0]+pad,plat[nb]+pad),
                      np.diff(plon)[-1]-2*pad, np.diff(plat)[-1]-2*pad,
                      fill=True, facecolor=inactive_color,
                      edgecolor='w'))
    ax.text(plon.mean(),nb, lab, fontsize=15,
             horizontalalignment='center', verticalalignment='center',
             color=tcol)

# label the buttons (numbered bottom to top, 0 to NB-1)
bdict = dict(zip(range(NB), blist))
Bdict = dict(zip(range(NB), Blist))
active_color = 'k'
inactive_color = 'w'
for bnum in bdict.keys():
    if bdict[bnum] == 'start':
        addButtonLabel(ax2, xbc, ybc, bnum, Bdict[bnum], tcol=active_color)
    else:
        addButtonLabel(ax2, xbc, ybc, bnum, Bdict[bnum], tcol=inactive_color)

plt.show()
pfun.topfig()

# polygon functions
def get_indices_in_polygon(plon_poly, plat_poly, NR, NC):
    # get indices of points inside a polygon
    V = np.ones((len(plon_poly),2))
    V[:,0] = plon_poly
    V[:,1] = plat_poly
    P = mpath.Path(V)
    # grid centers
    x = np.arange(NC)
    y = np.arange(NR)
    # matrix versions of grids
    X, Y = np.meshgrid(x,y)
    M, L = X.shape
    Rlon = X.flatten()
    Rlat = Y.flatten()
    R = np.ones((len(Rlon),2))
    R[:,0] = Rlon
    R[:,1] = Rlat
    RR = P.contains_points(R) # boolean
    # create arrays of i (column) and j (row) indices
    i_rho = np.arange(L).reshape((1,L)).repeat(M, axis=0)
    j_rho = np.arange(M).reshape((M,1)).repeat(L, axis=1)
    # pack indices that are inside the polygon
    # as a numpy int array, with two columns, packed in order j,i
    ji_rho_in = np.array([j_rho.flatten()[RR], i_rho.flatten()[RR]],
                         dtype=int).T
    return ji_rho_in
def remove_poly():    
    try: # remove old polygon lines if they exist
        pl = pline.pop(0)
        pl.remove()
    except (NameError, IndexError):
        pass

# allow user to edit mask until done
flag_get_ginput = True # Make False to exit the ginput loop
flag_continue = False # need to push START to make this True
flag_start = True # to ensure we only push the start button once
flag_e = 'p' # 'p' for polygon routines
pline = []
plon_poly = []
plat_poly = []

while flag_get_ginput:

    # get ginput, note that you can click with any key
    a = plt.ginput(n=1, timeout=-1)
    # returns a list of tuples - of length 1
    b = np.array(a)
    b = np.round(b).astype(int)
    if b.shape != (1,2):
        b = np.array([[-1000, -1000]])
    ix = b[0, 0]
    iy = b[0, 1]

    # this code deals with button input
    if (ix >= offset):
        # were are in the buttons
        nb = iy # button number
        if (bdict[nb]=='start') and flag_start:
            flag_start = False
            flag_continue = False
            cs.set_data(h)
            ax1.set_title('PAUSED: Initial Mask')
            # reset button colors
            for bnum in bdict.keys():
                if bdict[bnum] == 'start':
                    addButtonLabel(ax2, xbc, ybc, bnum, Bdict[bnum], tcol=inactive_color)
                else:
                    addButtonLabel(ax2, xbc, ybc, bnum, Bdict[bnum], tcol=active_color)
        elif (bdict[nb]=='pause') and not flag_start:
            flag_continue = False
            ax1.set_title('PAUSED')
        elif (bdict[nb]=='startPoly') and not flag_start:
            flag_continue = True
            flag_e = 'p'
            ax1.set_title('Click to Add Poly Points')
            remove_poly()
            pline = []
            plon_poly = []
            plat_poly = []
        elif (bdict[nb]=='polyInfo') and not flag_start:
            flag_continue = False
            ji_rho_in = get_indices_in_polygon(plon_poly, plat_poly, NR, NC)
            dap = da[ji_rho_in[:,0], ji_rho_in[:,1]]
            dvp = h[ji_rho_in[:,0], ji_rho_in[:,1]] * da[ji_rho_in[:,0], ji_rho_in[:,1]]
            ap = np.nansum(dap)
            vp = np.nansum(dvp)
            hp = vp/ap
            print('Volume inside polygon = %0.1f km3' % (vp/1e9) )
            print('Area inside polygon = %0.1f km2' % (ap/1e6) )
            print('Mean Depth inside polygon = %0.1f m' % (hp) )
            inp = input('Push Return to continue\n')
            remove_poly()
            ax1.set_title('PAUSED')
        elif (bdict[nb]=='lineLength') and not flag_start:
            flag_continue = False
            x = plon_poly
            y = plat_poly
            dist = 0
            for ii in range(len(x)-1):
                dist += np.sqrt( (xm[y[ii+1],x[ii+1]]-xm[y[ii],x[ii]])**2
                    + (ym[y[ii+1],x[ii+1]]-ym[y[ii],x[ii]])**2 )
            print('Line length = %0.1f km' % (dist/1e3))
            inp = input('Push Return to continue\n')
            remove_poly()
            ax1.set_title('PAUSED')
        elif (bdict[nb]=='done') and not flag_start:
            flag_get_ginput = False
            ax1.set_title('DONE', fontweight='bold', color='b')
        else:
            pass
        plt.draw()

    # this code deals with map input, and only responds when
    # we are clicking on the map
    elif flag_continue and not flag_start:
        # we are in the data field
        if flag_e == 'p':
            # this draws a polygon as you click
            plon_poly.append(ix)
            plat_poly.append(iy)
            remove_poly()
            aa = ax1.axis()
            pline = ax1.plot(plon_poly, plat_poly,'*-r')
            ax1.axis(aa)
            ax1.set_title('Adding to Polygon')
        plt.draw()

