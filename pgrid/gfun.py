# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 16:18:09 2016

@author: PM5

Utility function for pgrid.
"""

# USER EDIT

gridname = 'cas1'

dir0 = '/Users/PM5/Documents/'
pgdir = dir0 + 'ptools_output/pgrid/'

if 'aestus' in gridname:
    ri_dir = dir0 + 'ptools_output/river/analytical/'
else:
    ri_dir = dir0 + 'ptools_output/river/pnw_all_2016_07/'

# END USER EDIT

import os
import sys
alp = os.path.abspath(dir0 +'LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

plp = os.path.abspath(dir0 +'LiveOcean/plotting')
if plp not in sys.path:
    sys.path.append(plp)

import numpy as np
import pandas as pd

def gstart():
    gdir = pgdir + gridname + '/'
    Gr ={'gridname': gridname, 'dir0': dir0, 'pgdir': pgdir, 'gdir': gdir,
         'ri_dir': ri_dir}
    return Gr

def select_file(Gr, using_old_grid=False):
    # interactive selection
    if using_old_grid==True:
        fn_list = []
        dir0 = '/Users/PM5/Documents/LiveOcean_data/grids/'
        gn_list = ['cascadia1', 'cascadia2']
        for gn in gn_list:
            fn_list.append(dir0 + gn + '/grid.nc')
    elif using_old_grid==False:
        print('\n** %s in <<%s>> **\n' % ('Choose file to edit', Gr['gridname']))
        fn_list_raw = os.listdir(Gr['gdir'])
        fn_list = []
        for item in fn_list_raw:
            if item[-3:] == '.nc':
                fn_list.append(item)
    Nfn = len(fn_list)
    fn_dict = dict(zip(range(Nfn), fn_list))
    for nfn in range(Nfn):
        print(str(nfn) + ': ' + fn_list[nfn])
    my_nfn = int(input('-- Input number -- '))
    fn = fn_dict[my_nfn]
    return fn

def increment_filename(fn, tag='_m'):
    # create the new file name
    gni = fn.find(tag)
    new_num = ('00' + str(int(fn[gni+2: gni+4]) + 1))[-2:]
    fn_new = fn.replace(fn[gni:gni+4], tag + new_num)

    return fn_new

def GRID_PlusMinusScheme_rx0(MSK, Hobs, rx0max, AreaMatrix):
    """
    This is a faster version of GRID_PlusMinusScheme_rx0_ORIGr, with about 15x
    speedup in the 100x200 test grid.  It is comparable to the Matlab version.

    ** The depth matrix Hobs MUST BE POSITIVE in non-masked cells **

    The results were nearly identical to those from GRID_PlusMinusScheme_rx0,
    but with some variation up to +/- 45 m in some grid cells.  I suspect that
    this is do to the fact that the order in which I flip the grid around is
    different than in the original.  Since I see no reason for this order
    to be important I will assume the difference is not important.
    
    With fjord_cliff_edges=True is deviates from its usual volume-conserving
    nature when it is next to a masked region, and instead adjusts the slope
    by preferentially deepening at the coast.  This does a much better job of
    preserving thalweg depth in channels like Hood Canal.
    
    """
    fjord_cliff_edges = True
    
    HH=Hobs.copy()
    AA = AreaMatrix.copy()
    MM = MSK.copy()
    R=(1-rx0max)/(1+rx0max)
    tol=0.000001
    IsFinished = 1
    count = 0
    maxcount = 1000
    while True and count < maxcount:
        IsFinished=1
        for ff in range(5):
            if ff == 0:
                do_smooth = True
            elif ff == 1:
                do_smooth = True
                HH = np.fliplr(HH)
                AA = np.fliplr(AA)
                MM = np.fliplr(MM)
            elif ff == 2:
                do_smooth = True
                HH = HH.T
                AA = AA.T
                MM = MM.T
            elif ff == 3:
                do_smooth = True
                HH = np.fliplr(HH)
                AA = np.fliplr(AA)
                MM = np.fliplr(MM)
            elif ff == 4:
                do_smooth = False
                HH = HH.T
                HH = np.fliplr(HH)
                HH = np.flipud(HH)
                AA = AA.T
                AA = np.fliplr(AA)
                AA = np.flipud(AA)
                MM = MM.T
                MM = np.fliplr(MM)
                MM = np.flipud(MM)
            if do_smooth:
                NR, NC = HH.shape
                for ii in range(NC-1):
                    H = HH[:, ii]
                    Hn = HH[:, ii+1]
                    A = AA[:, ii]
                    An = AA[:, ii+1]
                    M = MM[:, ii]
                    Mn = MM[:, ii+1]
                    LowerBound = Hn*R
                    # mask is true when Hn is significantly deeper than H
                    # and when both are water points
                    # and when these are the case it makes H a little deeper
                    # and Hn a litte shallower
                    mask = (H - LowerBound < -tol) & (M == 1) & (Mn == 1)
                    if np.any(mask):
                        IsFinished=0
                        h = (R*Hn - H)/(An + R*A)
                        if ii > 0 and fjord_cliff_edges:
                            Mm = MM[:, ii-1]
                            xm = Mm == 0 # true when there is land to the left
                            xH = H.copy()
                            xH[xm] = xH[xm] + 2*An[xm]*h[xm]
                            xH[~xm] = xH[~xm] + An[~xm]*h[~xm]
                            H = xH.copy()
                            xHn = Hn.copy()
                            xHn[xm] = xHn[xm] - 0*A[xm]*h[xm]
                            xHn[~xm] = xHn[~xm] - A[~xm]*h[~xm]
                            Hn = xHn.copy()
                        else:
                            H = H + An*h
                            Hn = Hn - A*h
                        HH[mask, ii] = H[mask]
                        HH[mask, ii + 1] = Hn[mask]
                        
        if IsFinished == 1:
            break
        #print(count)
        count += 1
    print('Number of iterations = ' + str(count))
    if count == maxcount:
        print('\n** WARNING: more iterations needed! **\n')

    return HH
    
def get_plon_plat(using_old_grid, ds):
    if using_old_grid==True:
        # because older grids did not have lon,lat_psi_ex we create this
        # as an extension of lon,lat_psi
        plon0 = ds.variables['lon_psi'][:]
        plat0 = ds.variables['lat_psi'][:]
        dx = plon0[0,1] - plon0[0,0]
        dy = plat0[1,0] - plat0[0,0]
        ny0, nx0 = plon0.shape
        plon = np.nan * np.ones((ny0+2, nx0+2))
        plat = np.nan * np.ones((ny0+2, nx0+2))
        plon[1:-1, 1:-1] = plon0
        plat[1:-1, 1:-1] = plat0
        plon[:,0] = plon0[0,0] - dx
        plon[:,-1] = plon0[0,-1] + dx
        plon[0,:] = plon[1,:]
        plon[-1,:] = plon[-2,:]
        plat[0,:] = plat0[0,0] - dy
        plat[-1,:] = plat0[-1,0] + dy
        plat[:,0] = plat[:,1]
        plat[:,-1] = plat[:,-2]
    elif using_old_grid==False:
        plon = ds.variables['lon_psi_ex'][:]
        plat = ds.variables['lat_psi_ex'][:]
    return (plon, plat)
    
def get_grids(ds):
    lon_dict = dict()
    lat_dict = dict()
    mask_dict = dict()
    tag_list = ['rho', 'u', 'v', 'psi']
    for tag in tag_list:
        lon_dict[tag] = ds.variables['lon_'+tag][:]
        lat_dict[tag] = ds.variables['lat_'+tag][:]
        mask_dict[tag] = ds.variables['mask_'+tag][:]
    return (lon_dict, lat_dict, mask_dict)
    
def add_river_tracks(Gr, ds, ax):
    # add river tracks and endpoints
    in_rfn = Gr['gdir'] + 'river_info.csv'
    try:
        df = pd.read_csv(in_rfn, index_col='rname')
    except FileNotFoundError:
        return
    uv_dict = df['uv']
    row_dict_py = df['row_py']
    col_dict_py = df['col_py']
    isign_dict = df['isign']
    idir_dict = df['idir']
    # get grids
    lon_dict, lat_dict, mask_dict = get_grids(ds)
    lonu = lon_dict['u']
    latu = lat_dict['u']
    lonv = lon_dict['v']
    latv = lat_dict['v']
    # plot river tracks
    for rn in df.index:
        fn_tr = Gr['ri_dir'] + 'tracks/' + rn + '.csv'
        df_tr = pd.read_csv(fn_tr, index_col='ind')
        x = df_tr['lon'].values
        y = df_tr['lat'].values
        ax.plot(x, y, '-r', linewidth=2)
        ax.plot(x[-1], y[-1], '*r')    
        if uv_dict[rn] == 'u' and isign_dict[rn] == 1:
            ax.plot(lonu[row_dict_py[rn], col_dict_py[rn]],
                    latu[row_dict_py[rn], col_dict_py[rn]], '>r')
        elif uv_dict[rn] == 'u' and isign_dict[rn] == -1:
            ax.plot(lonu[row_dict_py[rn], col_dict_py[rn]],
                    latu[row_dict_py[rn], col_dict_py[rn]], '<r')
        elif uv_dict[rn] == 'v' and isign_dict[rn] == 1:
            ax.plot(lonv[row_dict_py[rn], col_dict_py[rn]],
                    latv[row_dict_py[rn], col_dict_py[rn]], '^b')
        elif uv_dict[rn] == 'v' and isign_dict[rn] == -1:
            ax.plot(lonv[row_dict_py[rn], col_dict_py[rn]],
                    latv[row_dict_py[rn], col_dict_py[rn]], 'vb')
                    
def edit_mask_river_tracks(Gr, NR, ax):
    # add river tracks and endpoints for edit_mask.py
    in_rfn = Gr['gdir'] + 'river_info.csv'
    try:
        df = pd.read_csv(in_rfn, index_col='rname')
    except FileNotFoundError:
        return
    uv_dict = df['uv']
    row_dict_py = df['row_py']
    col_dict_py = df['col_py']
    isign_dict = df['isign']
    idir_dict = df['idir']
    # plot river tracks
    for rn in df.index:
        yy = NR - row_dict_py[rn] - 1
        xx = col_dict_py[rn]
        if uv_dict[rn] == 'u' and isign_dict[rn] == 1:
            ax.plot(xx+.5, yy, '>r')
        elif uv_dict[rn] == 'u' and isign_dict[rn] == -1:
            ax.plot(xx+.5, yy, '<r')
        elif uv_dict[rn] == 'v' and isign_dict[rn] == 1:
            ax.plot(xx, yy-.5, '^b')
        elif uv_dict[rn] == 'v' and isign_dict[rn] == -1:
            ax.plot(xx, yy-.5, 'vb')    

def show_z_info(zm, ax):
    # find the max value of z (DEBUGGING)
    (rowmax, colmax) = np.unravel_index(np.argmax(zm), zm.shape)
    zmax = zm[rowmax, colmax]
    print('Max z = ' + str(zmax))
    lon_rho = ds['lon_rho'][:]
    lat_rho = ds['lat_rho'][:]
    ax.plot(lon_rho[rowmax, colmax], lat_rho[rowmax, colmax], '*m', markersize=20)

def show_grids(ds, ax):
    lon_dict, lat_dict, mask_dict = get_grids(ds)
    marker_dict = {'rho': 'ok',
                 'u': '>r',
                 'v': '^b',
                 'psi': 'xg'}
    for tag in tag_list:
        ax.plot(lon_dict[tag][mask_dict[tag]==1], lat_dict[tag][mask_dict[tag]==1],
                marker_dict[tag])
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(pfun.get_aa(ds))



