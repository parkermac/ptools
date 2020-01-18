"""
Code to test speed and scalability of 2-D interpolation.

Specifically I am benchmarking my own zfun.interp_scattered on plaid()
with a nearest neighbor algorithm.

In this case it is applied to the cas6 grid
"""

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
sys.path.append(os.path.abspath('../../LiveOcean/plotting'))
import pfun
import zfun
import zrfun
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
Ldir['gtagex'] = Ldir['gtag'] + '_lo8b'


import numpy as np
import numpy.random as rand
from scipy.spatial import cKDTree
from time import time
import pandas as pd
import matplotlib.pyplot as plt
import netCDF4 as nc

fn = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f2017.07.04/ocean_his_0001.nc'
G = zrfun.get_basic_info(fn, only_G=True)
ds = nc.Dataset(fn)
Z = ds['salt'][0,-1,:,:]

X = G['lon_rho']
Y = G['lat_rho']
xvec = X[0,:].flatten()
yvec = Y[:,0].flatten()

# create the nearest neighbor Tree object
# (the speed of the method relies on doing this)
Mask = ~Z.mask
xy = np.array((X[Mask],Y[Mask])).T
xyT = cKDTree(xy)
ZZ = Z[Mask].data

# the goal here is to find the value of Z at a bunch of arbitrary x,y locations
# and then compare different ways of making this interpolation

# ASSUMES we are working with a plaid grid

# list of the number of particle positions to loop over
S_list = [100, 1000, 10000, 100000]
# and a DataFrame to store the results (the interpolation times)
df = pd.DataFrame(index=S_list, columns=['interp2','nearest'])

for S in S_list:

    # create random particle positions to interpolate Z to
    xs = 2*rand.random_sample(S) - 125.5 # -125.5 for masked, -127 for unmasked
    ys = 2*rand.random_sample(S) + 46.5

    # (1) Use my zfun tools (bi-linear interpolation)
    tt0 = time()
    # note: zfun.interp_scattered_on_plaid() uses zfun.get_interpolant()
    zs = zfun.interp_scattered_on_plaid(xs, ys, xvec, yvec, Z)
    tt1 = time() - tt0

    # use nearest neighbor with the tree already created
    tt0 = time()
    xys = np.array((xs,ys)).T
    # use n_jobs=-1 to use all available cores
    zs_alt = ZZ[xyT.query(xys, n_jobs=-1)[1]]
    tt2 = time() - tt0

    df.loc[S,'interp2'] = tt1
    df.loc[S,'nearest'] = tt2

    print('S = %10d: interp2 = %0.4f, nearest = %0.4f' % (S, tt1, tt2))

plt.close('all')
fig = plt.figure(figsize=(13,5))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

df.plot(style='-o',ax=ax1)

ax2.plot(zs, zs_alt, '.c')
ax2.axis('square')
ax2.grid(True)

ax3.pcolormesh(X,Y,Z, vmin=zs_alt.min(), vmax=zs_alt.max(), cmap='rainbow')
ax3.plot(xs[::100], ys[::100], '.k', markersize=.3, alpha = .6)
pfun.dar(ax3)

plt.show()

