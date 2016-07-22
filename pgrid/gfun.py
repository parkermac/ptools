# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 16:18:09 2016

@author: PM5

Utility function for pgrid.
"""

# USER EDIT

gridname = 'cascadia2'
dir0 = '/Users/PM5/Documents/'
pgdir = dir0 + 'ptools_output/pgrid/'
ri_dir = dir0 + 'ptools_output/river/pnw_all_2016_07/'

# END USER EDIT

import os
import sys
# "alp" should also be on the path for any calling code
alp = os.path.abspath(dir0 +'LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)

import numpy as np
import matfun

def gstart():
    gdir = pgdir + gridname + '/'
    G = {'gridname': gridname, 'dir0': dir0, 'pgdir': pgdir, 'gdir': gdir,
         'ri_dir': ri_dir}
    return G

def select_file(G):
    # interactive selection
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

def get_coast():
    c_dir = dir0 + 'tools_data/geo_data/coast/'
    c_file = 'pnw_coast_combined.mat'
    c_fn = c_dir + c_file
    cmat = matfun.loadmat(c_fn)

    return cmat

def GRID_PlusMinusScheme_rx0(MSK, Hobs, rx0max, AreaMatrix):
    """
    This was coded by copying matlab code from LP_Bathymetry.
    It appears to give identical results to those from Matlab, but for
    some reason its performance is 10x slower.

    The depth matrix Hobs MUST BE POSITIVE outside of masked locations.

    """

    #% ---MSK(eta_rho,xi_rho) is the mask of the grid
    #%      1 for sea
    #%      0 for land
    #% ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
    #% ---rx0max is the target rx0 roughness factor
    #% ---AreaMatrix(eta_rho,xi_rho) is the matrix of areas at
    #% rho-points.

    eta_rho, xi_rho = Hobs.shape
    ListNeigh=np.array([[1, 0], [0, 1], [-1, 0], [0, -1]])
    RetBathy=Hobs.copy()

    HmodifVal=0
    TheMultiplier=(1-rx0max)/(1+rx0max)
    tol=0.000001
    ValueFct=0
    IsFinished = 1
    count = 0
    while True and count < 1000:
        if np.mod(count,10) == 0:
            print('Count = ' + str(count))
        count += 1
        IsFinished=1
        for iEta in range(eta_rho):
            for iXi in range(xi_rho):
                if MSK[iEta, iXi] == 1:
                    Area=AreaMatrix[iEta, iXi]
                    for ineigh in range(4):
                        iEtaN=iEta+ListNeigh[ineigh,0];
                        iXiN=iXi+ListNeigh[ineigh,1];
                        if ((iEtaN < eta_rho) and (iEtaN >= 0)
                            and (iXiN < xi_rho) and (iXiN >= 0)
                            and MSK[iEtaN, iXiN] == 1):
                            AreaN=AreaMatrix[iEtaN, iXiN]
                            LowerBound=RetBathy[iEtaN, iXiN]*TheMultiplier
                            if (RetBathy[iEta,iXi] - LowerBound < -tol):
                                IsFinished=0
                                h=(TheMultiplier*RetBathy[iEtaN, iXiN]-RetBathy[iEta, iXi])/(AreaN+TheMultiplier*Area)
                                RetBathy[iEta, iXi]=RetBathy[iEta, iXi]+AreaN*h
                                RetBathy[iEtaN, iXiN]=RetBathy[iEtaN, iXiN]-Area*h
                                HmodifVal=HmodifVal+np.abs(h)
                                ValueFct=ValueFct+np.abs(h)*(Area+AreaN)
        if IsFinished == 1:
            break

    H=AreaMatrix*Hobs*MSK
    TheBathymetry1=np.sum(H)
    H=AreaMatrix*RetBathy*MSK
    TheBathymetry2=np.sum(H)
    DeltaBathymetry=TheBathymetry1-TheBathymetry2
    print('Number of iterations = ' + str(count))
    print('DeltaBathymetry=' + str(DeltaBathymetry))

    return RetBathy, HmodifVal, ValueFct

def GRID_PlusMinusScheme_rx0_v2(MSK, Hobs, rx0max, AreaMatrix):
    """
    This is a faster version of GRID_PlusMinusScheme_rx0, with about 15x
    speedup in the 100x200 test grid.  It is comparable to the Matlab version.

    The results were nearly identical to those from GRID_PlusMinusScheme_rx0,
    but with some variation up to +/- 45 m in some grid cells.  I suspect that
    this is do to the fact that the order in which I flip the gird around is
    different than in the original.  Since I see no reason for this order
    to be important I will assume the difference is not important.

    The depth matrix Hobs MUST BE POSITIVE outside of masked locations.

    """

    #% ---MSK(eta_rho,xi_rho) is the mask of the grid
    #%      1 for sea
    #%      0 for land
    #% ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
    #% ---rx0max is the target rx0 roughness factor
    #% ---AreaMatrix(eta_rho,xi_rho) is the matrix of areas at
    #% rho-points.

    HH=Hobs.copy()

    AA = AreaMatrix.copy()

    MM = MSK.copy()

    R=(1-rx0max)/(1+rx0max)
    tol=0.000001
    IsFinished = 1
    count = 0
    maxcount = 1000
    while True and count < maxcount:
#        if np.mod(count,10) == 0:
#            print('Count = ' + str(count))
#        count += 1
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


