# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 16:45:54 2016

@author: PM5

Smooth the grid.

"""

from importlib import reload
import gfun; reload(gfun)
import pfun

import netCDF4 as nc
import os
import shutil
import numpy as np

gname = 'salish1/grid_m02_r01_s00_x00.nc'
in_fn = '/Users/PM5/Documents/ptools_output/pgrid/' + gname

# Test: retrieve the output

ds = nc.Dataset(in_fn)
z = -ds.variables['h'][:]
mask_rho = ds.variables['mask_rho'][:]
plon = ds.variables['lon_psi_ex'][:]
plat = ds.variables['lat_psi_ex'][:]
lon = ds.variables['lon_rho'][:]
lat = ds.variables['lat_rho'][:]
dx = 1/ds.variables['pm'][:]
dy = 1/ds.variables['pn'][:]

# Smoothing set up

# the land mask
# land = 0.
# water = 1.
MSK = mask_rho.copy()

# make sure that anything not masked is
# not shallower than min_depth
min_depth = 0. # a negative number would be like on land
Hobs_1 = -z.copy()
Hobs = Hobs_1.copy()
Hobs_1[Hobs_1 < min_depth] = min_depth
Hobs[mask_rho==1] = Hobs_1[mask_rho==1]
# and make it positive definite under water
offset = 1
Hobs = Hobs + offset
# offset is removed later

# smoothing target (rx0max <= 0.2 recommended)
rx0max = 0.15

# create the area matrix
AreaMatrix = dx * dy

#% create smoothed bathymetry

import time
tt0 = time.time()

#def GRID_PlusMinusScheme_rx0(MSK, Hobs, rx0max, AreaMatrix):

"""
This is a faster version of GRID_PlusMinusScheme_rx0_ORIG, with about 15x
speedup in the 100x200 test grid.  It is comparable to the Matlab version.

** The depth matrix Hobs MUST BE POSITIVE outside of masked locations **

The results were nearly identical to those from GRID_PlusMinusScheme_rx0,
but with some variation up to +/- 45 m in some grid cells.  I suspect that
this is do to the fact that the order in which I flip the grid around is
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
                # mask is true when Hn is significantly deeper than H
                # and when both are water points
                # and when these are the case it makes H a little deeper
                # and Hn a litte shallower
                mask = (H - LowerBound < -tol) & (M == 1) & (Mn == 1)
                if np.any(mask):
                    IsFinished=0
                    h = (R*Hn - H)/(An + R*A)
                    fjord_cliff_edges = True
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
        
    count += 1
    
print('Number of iterations = ' + str(count))
if count == maxcount:
    print('\n** WARNING: more iterations needed! **\n')

# correct for offset
HH = HH - offset

#return HH
Hnew = HH

#Hnew = gfun.GRID_PlusMinusScheme_rx0(MSK, Hobs, rx0max, AreaMatrix)

print('Smoothing took %0.1f seconds' % (time.time() - tt0))
zn = -Hnew

#%% plotting

if True:

    import matplotlib.pyplot as plt

    ax_lims = (plon[0,0], plon[0,-1], plat[0,0], plat[-1,0])

    z_m = np.ma.masked_where(mask_rho==0, z)
    zn_m = np.ma.masked_where(mask_rho==0, zn)
    dz = zn - z

    #plt.close()

    ds = nc.Dataset(in_fn)

    fig = plt.figure(figsize=(14,10))

    ax = fig.add_subplot(121)
    cmap1 = plt.get_cmap(name='terrain')
    cs = ax.pcolormesh(plon, plat, zn_m,
                       vmin=-200, vmax=200, cmap = cmap1)
    fig.colorbar(cs, ax=ax, extend='both')
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(pfun.get_aa(ds))
    ax.set_title(gname + ' NEW')

    ax = fig.add_subplot(122)
    cmap1 = plt.get_cmap(name='bwr')
    cs = ax.pcolormesh(plon, plat, dz,
                       vmin=-100, vmax=100, cmap = cmap1)
    fig.colorbar(cs, ax=ax, extend='both')
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(pfun.get_aa(ds))
    ax.set_title('New Z - Old Z (m)')

    plt.show()

    ds.close()
