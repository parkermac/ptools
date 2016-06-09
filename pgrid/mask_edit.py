# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 17:14:12 2016

@author: PM5

Code to edit a grid mask by hand.  We don't use masked arrays, but instead
keep one matrix Z that is the full bathymetry, and another M, which is
a boolean of the same size, True = masked.

In the editing we actually manipulate a flattened version "mi" that consists
of integer ones and zeros (1 for masked).  This is used to access an array
of the facecolors of all the tiles in the pcolormesh rendering of the
topography (NR*NC rows and 4 colums: each row is an RGBA color specification,
where A is transparency 0-1, smaller = more transparent).

When you edit the facecolor array it immediately updates the figure, and this
makes the operation relatively fast.  This also avoids the use of ax1.clear().
I tried an alternate version converting every thing to i,j indices and
then using imshow().  It was fast for just the topography, but bogged down
considerably when including the coast, which is essential.

The performance with a 300x600 grid gets a little slow, but is still OK.

"""

dir0 = '/Users/PM5/Documents/'

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as pltc
import shutil
import os
import sys
alp = os.path.abspath(dir0 + 'LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
from importlib import reload
import zfun; reload(zfun)
import matfun

flag_testing = False

if not flag_testing:

    # load bathymetry grid file

    indir = dir0 + 'ptools_output/pgrid/'

    if True:
        # interactive selection
        print('\n%s\n' % '** Choose file to edit **')
        fn_list_raw = os.listdir(indir)
        fn_list = []
        for item in fn_list_raw:
            if item[-3:] == '.nc':
                fn_list.append(item)
        Nfn = len(fn_list)
        fn_dict = dict(zip(range(Nfn), fn_list))
        for nfn in range(Nfn):
            print(str(nfn) + ': ' + fn_list[nfn])
        my_nfn = int(input('-- Input number -- '))
        gname = fn_dict[my_nfn]
    else:
        # test case created by make_grid.py
        gname = 'test_m00_r00_s00_x00.nc'

    in_fn = indir + gname

    # create the new file name
    gni = gname.find('_m')
    new_num = ('00' + str(int(gname[gni+2: gni+4]) + 1))[-2:]
    gname_new = gname.replace(gname[gni:gni+4],'_m' + new_num)
    out_fn = indir + gname_new

    # get the grid from NetCDF
    import netCDF4 as nc
    ds = nc.Dataset(in_fn)
    XC = ds.variables['lon_psi_ex'][:]
    YC = ds.variables['lat_psi_ex'][:]
    Z = ds.variables['z'][:]
    M = ds.variables['mask_rho'][:] == 0
    ds.close()
    xc = XC[0,:]
    yc = YC[:,0]
    # coastline
    do_coast = True
    if do_coast:
        c_dir = dir0 + 'tools_data/geo_data/coast/'
        c_file = 'pnw_coast_combined.mat'
        c_fn = c_dir + c_file
        cmat = matfun.loadmat(c_fn)

elif flag_testing:
    # simple grid for testing
    # grid corners
    xc = np.linspace(0,10,9)
    yc = np.linspace(0,10,11)
    # grid centers
    x = xc[:-1] + np.diff(xc)/2
    y = yc[:-1] + np.diff(yc)/2
    # matrix versions of grids
    X, Y = np.meshgrid(x,y)
    XC, YC = np.meshgrid(xc,yc)
    # data
    Z = (X**2 + Y**2) - 50
    # mask
    M = Z > 0 # True over land (positive Z)
    # coastline
    do_coast = False

# mi a vector version of the mask (True = masked) that we
# use when accessing the mesh colors
mi = M.flatten()
# this flattens by default in row-major (C-style) order, meaning it
# gives a row, appended by the next row, etc.

# size of the topo array
[NR, NC] = Z.shape

# PLOTTING
plt.close()

# set up the axes
fig = plt.figure(figsize=(10,10))
ax1 = plt.subplot2grid((1,3), (0,0), colspan=2)
ax2 = plt.subplot2grid((1,3), (0,2), colspan=1)

# initialize the data plot
cmap1 = plt.get_cmap(name='terrain')
cs = ax1.pcolormesh(XC,YC,Z, vmin=-200, vmax=200, cmap = cmap1)
if do_coast:
    ax1.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
    zfun.dar(ax1)
xl0 = (xc.min(), xc.max())
yl0 = (yc.min(), yc.max())
ax1.set_xlim(xl0)
ax1.set_ylim(yl0)
fig.colorbar(cs, ax=ax1, extend='both')

# create control buttons
NB = 4 # number of buttons
ybc = np.arange(NB+1)
offset = 1e5 # kludgey way to distinguish buttons from topography
xbc = np.arange(2) + offset
XBC, YBC = np.meshgrid(xbc, ybc)
ZB = np.arange(1,NB+1).reshape(NB,1)
cmap2 = plt.get_cmap(name='viridis')
ax2.pcolormesh(XBC,YBC,ZB, vmin=0, vmax=NB, cmap = cmap2)
ax2.set_axis_off()

def addButtonLabel(ax, xc, yc, nb, lab, tcol='k'):
    # draw and label buttons
    pad = .1
    ax.add_patch(
        plt.Rectangle((xc[0]+pad,yc[nb-1]+pad),
                      np.diff(xc)[-1]-2*pad, np.diff(yc)[-1]-2*pad,
                      fill=True, facecolor=inactive_color,
                      edgecolor='w'))
    ax.text(xc.mean(),yc[nb-1:nb+1].mean(), lab, fontsize=15,
             horizontalalignment='center', verticalalignment='center',
             color=tcol)

# label the buttons
NBstart = 1
NBpause = 2
NBcontinue = 3
NBdone = NB
active_color = 'k'
inactive_color = 'w'
mask_color = pltc.colorConverter.to_rgb('lavender')
addButtonLabel(ax2, xbc, ybc, NBpause, 'PAUSE', tcol=inactive_color)
addButtonLabel(ax2, xbc, ybc, NBcontinue, 'CONTINUE', tcol=inactive_color)
addButtonLabel(ax2, xbc, ybc, NBstart, 'START', tcol=active_color)
addButtonLabel(ax2, xbc, ybc, NBdone, 'DONE', tcol=inactive_color)

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
            col0 = ax1.collections[0]
            fc = col0.get_facecolor()
            fc0 = fc.copy() # store the original color
            for ii in range(len(mi)):
                # set flagged cells to mask_color
                if mi[ii] == True:
                    fc[ii,:3] = mask_color
            ax1.set_title('Initial Mask')
            # label the rest of the buttons
            addButtonLabel(ax2, xbc, ybc, NBstart, 'START', tcol=inactive_color)
            addButtonLabel(ax2, xbc, ybc, NBpause, 'PAUSE', tcol=active_color)
            addButtonLabel(ax2, xbc, ybc, NBcontinue, 'CONTINUE', tcol=active_color)
            addButtonLabel(ax2, xbc, ybc, NBdone, 'DONE', tcol=active_color)
            plt.draw()
        elif (np.ceil(b[:,1]).astype(int) == NBpause) and not flag_start:
            flag_continue = False
            ax1.set_title('PAUSED')
            plt.draw()
        elif (np.ceil(b[:,1]).astype(int) == NBcontinue) and not flag_start:
            flag_continue = True
            ax1.set_title('RESTARTED')
            plt.draw()
        elif (np.ceil(b[:,1]).astype(int) == NBdone) and not flag_start:
            flag_get_ginput = False
            ax1.set_xlim(xl0)
            ax1.set_ylim(yl0)
            ax1.set_title('DONE')
            plt.draw()
            # now copy the .nc file to a new name and then replace its mask.
        else:
            pass

    elif flag_continue and not flag_start:
        # we are in the data field
        ix0, ix1, frx = zfun.get_interpolant(b[:,0],xc)
        iy0, iy1, fry = zfun.get_interpolant(b[:,1],yc)
        iix = ix0[0]
        iiy = iy0[0]
        this_z = Z[iiy, iix]
        # this toggles the colors
        jj = iiy*NC + iix # index into the flattened array
        if mi[jj] == True:
            mi[jj] = False
            fc[jj,:3] = fc0[jj,:3] # original color
        elif mi[jj] == False:
            mi[jj] = True
            fc[jj,:3] = mask_color # 1 = white, 0 = black
        ax1.set_title('EDITING: iix=' + str(iix) + ' iiy=' + str(iiy)
                      + ' z=' + str(int(this_z)) + ' m')
        plt.draw()

#%% Save the output file

if not flag_testing:

    # get the original mask
    ds = nc.Dataset(in_fn)
    mask_rho_orig = ds['mask_rho'][:]
    ds.close()

    # create the new mask
    mii = mi.reshape(NR, NC)
    mask_rho = np.ones((NR,NC))
    mask_rho[mii == True] = 0.

    if not np.all(mask_rho == mask_rho_orig):
        print('Creating ' + gname_new)
        try:
            os.remove(out_fn)
        except OSError:
            pass # assume error was because the file did not exist
        shutil.copyfile(in_fn, out_fn)
        ds = nc.Dataset(out_fn, 'a')
        ds['mask_rho'][:] = mask_rho
        ds.close()
    else:
        print('No change to mask')
