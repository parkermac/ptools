"""
Code to test speed and scalability of 3-D interpolation using
nearest neighbor.  In this case it is applied to the cas6 grid.

This code pre-makes the KDTree(s) and saves them for use by
test_3d_interp_cas6.py.
"""

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
Ldir['gtagex'] = Ldir['gtag'] + '_lo8b'

import zrfun
import numpy.random as rand
from scipy.spatial import cKDTree
from time import time
import pickle
import netCDF4 as nc4
import numpy as np

outdir = Ldir['parent'] + 'ptools_output/tests/'
Lfun.make_dir(outdir)

# prepare fields to make the tree
tt0 = time()
fn = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f2017.07.04/ocean_his_0001.nc'
G, S, T = zrfun.get_basic_info(fn)
h = G['h']
zr = zrfun.get_z(h, 0*h, S, only_rho=True)
N,M,L = zr.shape
x = G['lon_rho']
y = G['lat_rho']
X = np.tile(x.reshape(1,M,L),[N,1,1])
Y = np.tile(y.reshape(1,M,L),[N,1,1])
H = np.tile(h.reshape(1,M,L),[N,1,1])
Z = zr/H # fractional depth (-1 to 0)
Mask = np.tile(G['mask_rho'].reshape(1,M,L),[N,1,1])
xyz = np.array((X[Mask],Y[Mask],Z[Mask])).T
print('Prepare fields to make tree %0.2f sec' % (time()-tt0))

# create the nearest neighbor Tree object
# (the speed of the method relies on doing this)
tt0 = time()
xyzT = cKDTree(xyz)
print('Create tree %0.2f sec' % (time()-tt0))

# save to disk
pickle.dump(xyzT, open(outdir + 'cas6_xyzT.p', 'wb'))