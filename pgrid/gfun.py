# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 16:18:09 2016

@author: PM5

Utility function for pgrid.
"""

# USER EDIT

gridname = 'bigSalish1'

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

def gstart():
    gdir = pgdir + gridname + '/'
    G = {'gridname': gridname, 'dir0': dir0, 'pgdir': pgdir, 'gdir': gdir,
         'ri_dir': ri_dir}
    return G

def select_file(G, LO=False):
    # interactive selection
    if LO==True:
        # use this to look at the final grids we have created
        fn_list = []
        dir0 = '/Users/PM5/Documents/LiveOcean/preamble/make_resources/'
        gn_list = ['cascadia1', 'cascadia2']
        for gn in gn_list:
            fn_list.append(dir0 + gn + '/grid.nc')
    elif LO==False:
        # use this to look at grids from pgrid
        print('\n** %s in <<%s>> **\n' % ('Choose file to edit', G['gridname']))
        fn_list_raw = os.listdir(G['gdir'])
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
    This is a faster version of GRID_PlusMinusScheme_rx0_ORIG, with about 15x
    speedup in the 100x200 test grid.  It is comparable to the Matlab version.

    ** The depth matrix Hobs MUST BE POSITIVE outside of masked locations **

    The results were nearly identical to those from GRID_PlusMinusScheme_rx0,
    but with some variation up to +/- 45 m in some grid cells.  I suspect that
    this is do to the fact that the order in which I flip the gird around is
    different than in the original.  Since I see no reason for this order
    to be important I will assume the difference is not important.
    """
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
                    mask = (H - LowerBound < -tol) & (M == 1) & (Mn == 1)
                    if np.any(mask):
                        IsFinished=0
                        h = (R*Hn - H)/(An + R*A)
                        H = H + An*h
                        Hn = Hn - A*h
                        HH[mask, ii] = H[mask]
                        HH[mask, ii + 1] = Hn[mask]
        if IsFinished == 1:
            break
    print('Number of iterations = ' + str(count))
    if count == maxcount:
        print('\n** WARNING: more iterations needed! **\n')

    return HH



