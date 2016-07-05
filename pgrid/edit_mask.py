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

from importlib import reload
import gfun; reload(gfun)
G = gfun.gstart()

import numpy as np
import shutil
import os
import pandas as pd

import zfun

import matplotlib.pyplot as plt
import matplotlib.colors as pltc
import sys
alp = os.path.abspath(G['dir0'] + 'LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
from importlib import reload

flag_testing = False

if not flag_testing:

    # select grid file
    fn = gfun.select_file(G)
    in_fn = G['gdir'] + fn
    # create new file name
    fn_new = gfun.increment_filename(fn, tag='_m')
    out_fn = G['gdir'] + fn_new

    # get the grid from NetCDF
    import netCDF4 as nc
    ds = nc.Dataset(in_fn)
    plon = ds.variables['lon_psi_ex'][:]
    plat = ds.variables['lat_psi_ex'][:]
    lonu = ds.variables['lon_u'][:]
    latu = ds.variables['lat_u'][:]
    lonv = ds.variables['lon_v'][:]
    latv = ds.variables['lat_v'][:]
    z = -ds.variables['h'][:]
    m = ds.variables['mask_rho'][:] == 0
    ds.close()
    plon = plon[0,:]
    plat = plat[:,0]
    # coastline
    do_coast = True
    if do_coast:
        cmat = gfun.get_coast()

elif flag_testing:
    # simple grid for testing
    # grid corners
    plon = np.linspace(0,10,9)
    plat = np.linspace(0,10,11)
    # grid centers
    x = plon[:-1] + np.diff(plon)/2
    y = plat[:-1] + np.diff(plat)/2
    # matrix versions of grids
    X, Y = np.meshgrid(x,y)
    plon, plat = np.meshgrid(plon,plat)
    # data
    z = (X**2 + Y**2) - 50
    # mask
    m = z > 0 # True over land (positive Z)
    # coastline
    do_coast = False

# mi a vector version of the mask (True = masked) that we
# use when accessing the mesh colors
mi = m.flatten()
# this flattens by default in row-major (C-style) order, meaning it
# gives a row, appended by the next row, etc.

# size of the topo array
[NR, NC] = z.shape

# PLOTTING
plt.close()

# set up the axes
fig = plt.figure(figsize=(10,10))
ax1 = plt.subplot2grid((1,3), (0,0), colspan=2)
ax2 = plt.subplot2grid((1,3), (0,2), colspan=1)

# initialize the data plot
cmap1 = plt.get_cmap(name='terrain')
tvmin = -200
tvmax = 200
cs = ax1.pcolormesh(plon,plat,z, vmin=tvmin, vmax=tvmax, cmap = cmap1)
if do_coast:
    ax1.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
    zfun.dar(ax1)

# add rivers
if not flag_testing:
    in_rfn = G['gdir'] + 'river_info.p'
    df = pd.read_pickle(in_rfn)

    for rn in df.index:
        fn_tr = G['ri_dir'] + 'tracks/' + rn + '.csv'
        df_tr = pd.read_csv(fn_tr, index_col='ind')
        x = df_tr['lon'].values
        y = df_tr['lat'].values
        ax1.plot(x, y, '-r', linewidth=2)
        ax1.plot(x[-1], y[-1], '*r')

        if df.ix[rn, 'uv'] == 'u' and df.ix[rn, 'isign'] == 1:
            ax1.plot(lonu[df.ix[rn, 'row_py'], df.ix[rn, 'col_py']],
                    latu[df.ix[rn, 'row_py'], df.ix[rn, 'col_py']], '>r')
        elif df.ix[rn, 'uv'] == 'u' and df.ix[rn, 'isign'] == -1:
            ax1.plot(lonu[df.ix[rn, 'row_py'], df.ix[rn, 'col_py']],
                    latu[df.ix[rn, 'row_py'], df.ix[rn, 'col_py']], '<r')
        elif df.ix[rn, 'uv'] == 'v' and df.ix[rn, 'isign'] == 1:
            ax1.plot(lonv[df.ix[rn, 'row_py'], df.ix[rn, 'col_py']],
                    latv[df.ix[rn, 'row_py'], df.ix[rn, 'col_py']], '^b')
        elif df.ix[rn, 'uv'] == 'v' and df.ix[rn, 'isign'] == -1:
            ax1.plot(lonv[df.ix[rn, 'row_py'], df.ix[rn, 'col_py']],
                    latv[df.ix[rn, 'row_py'], df.ix[rn, 'col_py']], 'vb')

xl0 = (plon.min(), plon.max())
yl0 = (plat.min(), plat.max())
ax1.set_xlim(xl0)
ax1.set_ylim(yl0)

fig.colorbar(cs, ax=ax1, extend='both')

# create control buttons
NB = 5 # number of buttons
ybc = np.arange(NB+1)
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
        plt.Rectangle((plon[0]+pad,plat[nb-1]+pad),
                      np.diff(plon)[-1]-2*pad, np.diff(plat)[-1]-2*pad,
                      fill=True, facecolor=inactive_color,
                      edgecolor='w'))
    ax.text(plon.mean(),plat[nb-1:nb+1].mean(), lab, fontsize=15,
             horizontalalignment='center', verticalalignment='center',
             color=tcol)

# label the buttons
NBstart = 1
NBpause = 2
NBcontinueM = 3
NBcontinueZ = 4
NBdone = NB
active_color = 'k'
inactive_color = 'w'
mask_color = pltc.colorConverter.to_rgb('lightsalmon')
addButtonLabel(ax2, xbc, ybc, NBpause, 'PAUSE', tcol=inactive_color)
addButtonLabel(ax2, xbc, ybc, NBcontinueM, 'CONTINUE\nEdit Mask', tcol=inactive_color)
addButtonLabel(ax2, xbc, ybc, NBcontinueZ, 'CONTINUE\nEdit Z', tcol=inactive_color)
addButtonLabel(ax2, xbc, ybc, NBstart, 'START', tcol=active_color)
addButtonLabel(ax2, xbc, ybc, NBdone, 'DONE', tcol=inactive_color)

plt.show()

# allow user to edit mask until done
flag_get_ginput = True # Make False to exit the ginput loop
flag_continue = False # need to push START to make this True
flag_start = True # to ensure we only push the start button once
flag_mz = 'm' # 'm' to edit mask, 'z' to edit z

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
            addButtonLabel(ax2, xbc, ybc, NBcontinueM, 'CONTINUE\nEdit Mask', tcol=active_color)
            addButtonLabel(ax2, xbc, ybc, NBcontinueZ, 'CONTINUE\nEdit Z', tcol=active_color)
            addButtonLabel(ax2, xbc, ybc, NBdone, 'DONE', tcol=active_color)
            plt.draw()
        elif (np.ceil(b[:,1]).astype(int) == NBpause) and not flag_start:
            flag_continue = False
            ax1.set_title('PAUSED')
            plt.draw()
        elif (np.ceil(b[:,1]).astype(int) == NBcontinueM) and not flag_start:
            flag_continue = True
            flag_mz = 'm'
            ax1.set_title('EDITING Mask')
            plt.draw()
        elif (np.ceil(b[:,1]).astype(int) == NBcontinueZ) and not flag_start:
            flag_continue = True
            flag_mz = 'z'
            ax1.set_title('EDITING Z')
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
        ix0, ix1, frx = zfun.get_interpolant(b[:,0],plon)
        iy0, iy1, fry = zfun.get_interpolant(b[:,1],plat)
        iix = ix0[0]
        iiy = iy0[0]
        this_z = z[iiy, iix]
        if flag_mz == 'm':
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
        elif flag_mz == 'z':
            # this carves to a specified depth, and removes the mask
            jj = iiy*NC + iix # index into the flattened array
            mi[jj] = False
            i_newcolor = np.arange(cmap1.N)
            z_newcolor = np.linspace(tvmin, tvmax, cmap1.N)
            ii_newcolor = i_newcolor[z_newcolor >= -5][0]
            fc_new = cmap1(ii_newcolor)[:3]
            fc[jj,:3] = fc_new # new color
            z[iiy, iix] = -5
            this_z = z[iiy, iix]
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
    m = mi.reshape(NR, NC)
    mask_rho = np.ones((NR,NC))
    mask_rho[m == True] = 0.

    if not np.all(mask_rho == mask_rho_orig):
        print('Creating ' + out_fn)
        try:
            os.remove(out_fn)
        except OSError:
            pass # assume error was because the file did not exist
        shutil.copyfile(in_fn, out_fn)
        ds = nc.Dataset(out_fn, 'a')
        ds['mask_rho'][:] = mask_rho
        ds['h'][:] = -z
        ds.close()
    else:
        print('No change to mask')
