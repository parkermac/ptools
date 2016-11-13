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

import gfun
from importlib import reload
reload(gfun)
G = gfun.gstart()
# running gfun.gstart() sets the path to include pfun and zfun
import pfun
import zfun

import numpy as np
import shutil
import os
import pandas as pd
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import matplotlib.colors as pltc
from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

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

elif flag_testing:
    # simple grid for testing
    # grid corners (like the psi_ex grid)
    plon = np.linspace(0,10,9) # 9
    plat = np.linspace(0,10,11) # 11
    # grid centers
    x = plon[:-1] + np.diff(plon)/2
    y = plat[:-1] + np.diff(plat)/2
    # matrix versions of grids
    X, Y = np.meshgrid(x,y)
    # synthetic topo data
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

#%% PLOTTING
# set up the axes
plt.close()
fig = plt.figure(figsize=(14,14))
ax1 = plt.subplot2grid((1,3), (0,0), colspan=2) # map
ax2 = plt.subplot2grid((1,3), (0,2), colspan=1) # buttons

#%% initialize the data plot
cmap1 = plt.get_cmap(name='terrain')
tvmin = -200
tvmax = 200
cs = ax1.pcolormesh(plon,plat,z, vmin=tvmin, vmax=tvmax, cmap = cmap1)
if do_coast:
    pfun.add_coast(ax1)
    pfun.dar(ax1)
# add rivers
if not flag_testing:
    in_rfn = G['gdir'] + 'river_info.csv'
    df = pd.read_csv(in_rfn, index_col='rname')
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
# set limits and colorbar
map_lims = [plon.min(), plon.max(), plat.min(), plat.max()]
ax1.axis(map_lims)
fig.colorbar(cs, ax=ax1, extend='both')

#%% create control buttons

# list is organized from bottom to top

blist = ['start', 'pause', 'continueM', 'continueZ',
         'polyToLand', 'polyToWater', 'startPoly',
         'done']

NB = len(blist) # number of buttons
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

# label the buttons (numbered bottom to top, 1 to NB)

bdict = dict(zip(range(1,NB+1), blist))

active_color = 'k'
inactive_color = 'w'
mask_color = pltc.colorConverter.to_rgb('lightsalmon')

for bnum in bdict.keys():
    if bdict[bnum] == 'start':
        addButtonLabel(ax2, xbc, ybc, bnum, bdict[bnum], tcol=active_color)
    else:
        addButtonLabel(ax2, xbc, ybc, bnum, bdict[bnum], tcol=inactive_color)

plt.show()
pfun.topfig()

#%% polygon functions

def get_indices_in_polygon(plon_poly, plat_poly, plon, plat):
    # get indices of points inside a polygon
    V = np.ones((len(plon_poly),2))
    V[:,0] = plon_poly
    V[:,1] = plat_poly
    P = mpath.Path(V)

    # grid centers
    # (plon and plat are vectors)
    x = plon[:-1] + np.diff(plon)/2
    y = plat[:-1] + np.diff(plat)/2
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

#%% allow user to edit mask until done
flag_get_ginput = True # Make False to exit the ginput loop
flag_continue = False # need to push START to make this True
flag_start = True # to ensure we only push the start button once
flag_e = 'm' # 'm' to edit mask, 'z' to edit z, 'p' for polygon routines
pline = []
plon_poly = []
plat_poly = []

while flag_get_ginput:

    # get ginput, note that you can click with any key
    a = plt.ginput(n=1, timeout=-1)
    # returns a list of tuples - of length 1
    b = np.array(a)
    if b.shape != (1,2):
        b = np.array([[-1000, -1000]])

    # this code deals with button input
    if (b[0,0] >= offset):
        # were are in the buttons
        nb = np.ceil(b[:,1]).astype(int)[0] # button number
        
        if (bdict[nb]=='start') and flag_start:
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
            # reset button colors           
            for bnum in bdict.keys():
                if bdict[bnum] == 'start':
                    addButtonLabel(ax2, xbc, ybc, bnum, bdict[bnum], tcol=inactive_color)
                else:
                    addButtonLabel(ax2, xbc, ybc, bnum, bdict[bnum], tcol=active_color)            
        elif (bdict[nb]=='pause') and not flag_start:
            flag_continue = False
            ax1.set_title('PAUSED')
        elif (bdict[nb]=='continueM') and not flag_start:
            flag_continue = True
            flag_e = 'm'
            ax1.set_title('EDITING Mask')
        elif (bdict[nb]=='continueZ') and not flag_start:
            flag_continue = True
            flag_e = 'z'
            ax1.set_title('EDITING Z')
        elif (bdict[nb]=='startPoly') and not flag_start:
            flag_continue = True
            flag_e = 'p'
            ax1.set_title('Click to Add Poly Points')
            remove_poly()
            pline = []
            plon_poly = []
            plat_poly = []
        elif (bdict[nb]=='polyToLand') and not flag_start:
            flag_continue = False
            ax1.set_title('Changing Poly to Land')
            ji_rho_in = get_indices_in_polygon(plon_poly, plat_poly, plon, plat)
            jj = ji_rho_in[:,0]*NC + ji_rho_in[:,1]
            mi[jj] = True
            fc[jj,:3] = mask_color
            remove_poly()
        elif (bdict[nb]=='polyToWater') and not flag_start:
            flag_continue = False
            ax1.set_title('Changing Poly to Water')
            ji_rho_in = get_indices_in_polygon(plon_poly, plat_poly, plon, plat)
            jj = ji_rho_in[:,0]*NC + ji_rho_in[:,1]
            mi[jj] = False
            fc[jj,:3] = fc0[jj,:3]
            remove_poly()
        elif (bdict[nb]=='done') and not flag_start:
            flag_get_ginput = False
            ax1.axis(map_lims)
            ax1.set_title('DONE')
        else:
            pass
        plt.draw()

    # this code deals with map input, and only responds when
    # we are clicking on the map
    elif flag_continue and not flag_start:
        # we are in the data field
        ix0, ix1, frx = zfun.get_interpolant(np.array(b[0,0]),plon)
        iy0, iy1, fry = zfun.get_interpolant(np.array(b[0,1]),plat)
        iix = ix0[0]
        iiy = iy0[0]
        this_z = z[iiy, iix]
        if flag_e == 'm':
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
        elif flag_e == 'z':
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
        elif flag_e == 'p':
            # this draws a polygon as you click
            plon_poly.append(b[0,0])
            plat_poly.append(b[0,1])
            remove_poly()
            pline = ax1.plot(plon_poly, plat_poly,'*-r')
            ax1.set_title('Adding to Polygon')
                       
        plt.draw()
        # somewhere I saw that the draw_idle command might be better for
        # redrawing, but I haven't seen a real difference
        #fig.canvas.draw_idle()

#%% Save the output file: executes when you push done button

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
