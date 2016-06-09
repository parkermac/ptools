# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 16:23:54 2016

@author: PM5

Code to edit a masked grid by hand,
using imshow instead of pcolormesh.

RESULT: This was faster then the original mask_edit.py IF I did not use the
coastline.  With the coastline on it was actually slower for things like
resizing and zooming.

"""

dir0 = '/Users/PM5/Documents/'

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
alp = os.path.abspath(dir0 + 'LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
from importlib import reload
import zfun; reload(zfun)
import matfun

if False:
    # simple grid for testing
    # grid corners
    gnx = 8
    gny = 10
    xc = np.linspace(0,gnx,gnx+1) - .5
    yc = np.linspace(0,gny,gny+1) - .5
    # grid centers
    x = xc[:-1] + np.diff(xc)/2
    y = yc[:-1] + np.diff(yc)/2
    # matrix versions of grids
    X, Y = np.meshgrid(x,y)
    # data
    Z = (X**2 + Y**2) - 50
    # mask
    M = Z > 0 # True over land (positive Z)
    do_coast = False
else:
    # actual bathymetry grid file
    indir = dir0 + '/ptools_output/pgrid/'
    inname = 'topo_test'
    in_fn = indir + inname + '.nc'
    # get the grid from NetCDF
    import netCDF4 as nc
    ds = nc.Dataset(in_fn)
    Z = ds.variables['z'][:]
    M = ds.variables['mask_rho'][:] == 1

    XC = ds.variables['lon_psi'][:]
    YC = ds.variables['lat_psi'][:]
    xcd = XC[0,:]
    ycd = YC[:,0]

    ds.close()
    gny, gnx = Z.shape
    xc = np.linspace(0,gnx,gnx+1) - .5
    yc = np.linspace(0,gny,gny+1) - .5
    do_coast = True

Z0 = Z.copy()
Z[M] = np.nan

if do_coast:
    # coastline
    # need to translate it into i,j space
    c_dir = dir0 + 'tools_data/geo_data/coast/'
    c_file = 'pnw_coast_combined.mat'
    c_fn = c_dir + c_file
    cmat = matfun.loadmat(c_fn)
    clon = cmat['lon']
    clat = cmat['lat']
    Clon = clon * np.nan
    Clat = clat * np.nan
    cmask = ~np.isnan(clon)
    clon_good = clon[cmask]
    clat_good = clat[cmask]
    ix0, ix1, frx = zfun.get_interpolant(clon_good,xcd)
    iy0, iy1, fry = zfun.get_interpolant(clat_good,ycd)
    cix = xc[ix0]*(1-frx) + xc[ix1]*frx
    ciy = yc[iy0]*(1-fry) + yc[iy1]*fry
    ZC = Z0.copy() * np.nan
    ZC[iy0, ix0] = 1
    Clon[cmask] = cix
    Clat[cmask] = ciy

# size of the topo array
[NR, NC] = Z.shape

# PLOTTING
plt.close()

# set up the axes
fig = plt.figure(figsize=(8,10))
ax1 = plt.subplot2grid((1,3), (0,0), colspan=2)
ax2 = plt.subplot2grid((1,3), (0,2), colspan=1)

# initialize the data plot
cmap1 = plt.get_cmap(name='Set3')
zmin = np.min(Z.flatten())
ax1.imshow(Z, interpolation='nearest')
if do_coast:
    ax1.plot(Clon,Clat,'-k')
    #ax1.imshow(ZC, interpolation='nearest')
xl0 = (xc.min(), xc.max())
yl0 = (yc.min(), yc.max())
xl = xl0
yl = yl0
ax1.set_xlim(xl0)
ax1.set_ylim(yl0)

# create control buttons
NB = 4 # number of buttons
ybc = np.arange(NB+1)
offset = 1e5 # kludgey way to distinguish buttons from topography
xbc = np.arange(2) + offset
XBC, YBC = np.meshgrid(xbc, ybc)
ZB = np.arange(1,NB+1).reshape(NB,1)
cmap2 = plt.get_cmap(name='Set3')
ax2.pcolormesh(XBC,YBC,ZB, vmin=0, vmax=NB, cmap = cmap2)
ax2.set_axis_off()

def addButtonLabel(ax, xc, yc, nb, lab, tcolor='k'):
    # draw and label buttons
    pad = .1
    ax.add_patch(
        plt.Rectangle((xc[0]+pad,yc[nb-1]+pad),
                      np.diff(xc)[-1]-2*pad, np.diff(yc)[-1]-2*pad,
                      fill=True, facecolor=inactive_color,
                      edgecolor='w'))
    ax.text(xc.mean(),yc[nb-1:nb+1].mean(), lab, fontsize=15,
             horizontalalignment='center', verticalalignment='center',
             color=tcolor)

# label the buttons
NBstart = 1
NBpause = 2
NBcontinue = 3
NBdone = NB
active_color = 'k'
inactive_color = 'lightgray'
addButtonLabel(ax2, xbc, ybc, NBpause, 'PAUSE', tcolor=inactive_color)
addButtonLabel(ax2, xbc, ybc, NBcontinue, 'CONTINUE', tcolor=inactive_color)
addButtonLabel(ax2, xbc, ybc, NBstart, 'START', tcolor=active_color)
addButtonLabel(ax2, xbc, ybc, NBdone, 'DONE', tcolor=inactive_color)

plt.show()

# allow user to edit mask until done
flag_get_ginput = True # Make False to exit the ginput loop
flag_continue = False # need to push START to make this True
flag_start = True # to ensure we only push the start button once

while flag_get_ginput:

    # get ginput
    a = plt.ginput(n=1) # returns a list of tuples - of length 1
    b = np.array(a)

    if (b[:,0] >= offset):
        # were are in the buttons
        if (np.ceil(b[:,1]).astype(int) == NBstart) and flag_start:
            flag_start = False
            flag_continue = True
            flag_get_colors = False
            ax1.set_title('Initial Mask')
            addButtonLabel(ax2, xbc, ybc, NBpause, 'PAUSE', tcolor=active_color)
            addButtonLabel(ax2, xbc, ybc, NBcontinue, 'CONTINUE', tcolor=active_color)
            addButtonLabel(ax2, xbc, ybc, NBstart, 'START', tcolor=inactive_color)
            addButtonLabel(ax2, xbc, ybc, NBdone, 'DONE', tcolor=active_color)
            plt.draw()
        elif np.ceil(b[:,1]).astype(int) == NBpause:
            flag_continue = False
            ax1.set_title('PAUSED')
            plt.draw()
        elif np.ceil(b[:,1]).astype(int) == NBcontinue:
            flag_continue = True
            xl = ax1.get_xlim()
            yl = ax1.get_ylim()
            ax1.set_title('RESTARTED')
            plt.draw()
        elif np.ceil(b[:,1]).astype(int) == NBdone:
            flag_get_ginput = False
            ax1.set_xlim(xl0)
            ax1.set_ylim(yl0)
            ax1.set_title('DONE')
            plt.draw()
        else:
            pass

    elif flag_continue:
        # we are in the data field
        ix0, ix1, frx = zfun.get_interpolant(b[:,0],xc)
        iy0, iy1, fry = zfun.get_interpolant(b[:,1],yc)
        iix = ix0[0]
        iiy = iy0[0]
        # continue
        # this toggles the colors
        if np.isnan(Z[iiy,iix]):
#            print('found a nan')
#            sys.stdout.flush()
            Z[iiy,iix] = Z0[iiy,iix] # original color
        else:
#            print('found a good point')
#            sys.stdout.flush()
            Z[iiy,iix] = np.nan
        ax1.clear()
        ax1.imshow(Z, interpolation='nearest')
        if do_coast:
            ax1.plot(Clon,Clat,'-k')
            #ax1.imshow(ZC, interpolation='nearest')
        ax1.set_title('EDITING: iix=' + str(iix) + ' iiy=' + str(iiy))
        ax1.set_xlim(xl)
        ax1.set_ylim(yl)
        plt.draw()
