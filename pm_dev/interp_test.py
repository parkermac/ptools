"""
Code to test interpolation schemes.

"""

# setup
import os
import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()
from importlib import reload
import zfun
reload(zfun)

import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

# this will be a test of filling out a plaid grid from
# scattered data, using nearest neighbor interpolation

# the plaid grid of box edges
xx = np.linspace(-10,10,21)
yy= xx
XX, YY = np.meshgrid(xx, yy)

# the plaid grid of box centers
x = xx[:-1] + np.diff(xx)/2
y = x
X, Y = np.meshgrid(x, y)

# data function
def get_z(x, y):
    # inputs are numpy arrays of any shape
    z = np.exp(-(x/7)**2 - (y/7)**2) 
    return z

# create a perfect data field on box centers
Z = get_z(X, Y)

if True:
    # this makes a data field at box centers starting from the
    # perfect data field but then masking out some of it
    fld = np.ma.masked_where(X+Y > 0, Z)
    # do the extrapolation using nearest neighbor
    fldf = fld.copy() # initialize the "filled" field
    xyorig = np.array((X[~fld.mask],Y[~fld.mask])).T
    xynew = np.array((X[fld.mask],Y[fld.mask])).T
    a = cKDTree(xyorig).query(xynew)
    aa = a[1]
    fldf[fld.mask] = fld[~fld.mask][aa]
    plot_points = False
else:
    # this makes scattered data
    fld = Z
    npts = 200
    xs = 20*np.random.rand(npts) - 10
    ys = 20*np.random.rand(npts) - 10
    zs = np.exp(-(xs/7)**2 - (ys/7)**2) 
    # do the extrapolation using nearest neighbor
    fldf = np.nan * np.ones(X.shape) # initialize the "filled" field
    xyorig = np.array((xs,ys)).T # data locations
    xynew = np.array((X,Y)).T # where I want to have the data
    a = cKDTree(xyorig).query(xynew)
    aa = a[1]
    # aa is an array the indices into the good data that is nearest
    # to each of the new locations
    fldf = zs[aa]  
    plot_points = True               

# PLOTTING

plt.close()

fig = plt.figure(figsize=(14,5))

ax = fig.add_subplot(121)
ax.set_aspect('equal')
cs = ax.pcolormesh(XX, YY, fld, vmin=0, vmax=1)
if plot_points:
    ax.plot(xs, ys, '.k')
fig.colorbar(cs)
ax.set_title('Original Field')

ax = fig.add_subplot(122)
ax.set_aspect('equal')
cs = ax.pcolormesh(XX, YY, fldf, vmin=0, vmax=1)
fig.colorbar(cs)
ax.set_title('Gridded/Extrapolated Field')

plt.show()