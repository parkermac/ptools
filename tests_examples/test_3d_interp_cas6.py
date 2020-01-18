"""
Code to test speed and scalability of 3-D interpolation using
nearest neighbor.  The goal here is to find the value of D at
a bunch of arbitrary x,y,z locations. In this case it is applied to the cas6 grid.

RESULT: the performance looks great to me
NP =        100: 0.0012 sec
NP =       1000: 0.0039 sec
NP =      10000: 0.0279 sec
NP =     100000: 0.2672 sec

"""

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
Ldir['gtagex'] = Ldir['gtag'] + '_lo8b'

import zrfun
import numpy as np
import numpy.random as rand
from time import time
import matplotlib.pyplot as plt
import netCDF4 as nc4
import pickle

# prepare data field to query
fn = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f2017.07.04/ocean_his_0001.nc'
if False:
    ds = nc4.Dataset(fn)
    D = ds['salt'][0,:,:,:] # the data field
    do_plot = False
else:
    # by making the data equal to the scaled depth we can
    # test the search strategy by comparing the output "dp"
    # to the input scaled depth "zp"
    G, S, T = zrfun.get_basic_info(fn)
    h = G['h']
    zr = zrfun.get_z(h, 0*h, S, only_rho=True)
    N,M,L = zr.shape
    H = np.tile(h.reshape(1,M,L),[N,1,1])
    Z = zr/H # scaled depth (-1 to 0)
    Mask = np.tile(G['mask_rho'].reshape(1,M,L),[N,1,1])
    D = Z.copy()
    do_plot = True
DD = D[Mask].data

# load the tree
outdir = Ldir['parent'] + 'ptools_output/tests/'
xyzT = pickle.load(open(outdir + 'cas6_xyzT.p', 'rb'))

# list of the number of particle positions to loop over
NP_list = [100, 1000, 10000]

for NP in NP_list:

    # create random particle positions to interpolate Z to
    xp = 2*rand.random_sample(NP) - 125.5 # -125.5 for masked, -127 for unmasked
    yp = 2*rand.random_sample(NP) + 46.5
    zp = rand.random_sample(NP) - 1
    xyzs = np.array((xp,yp,zp)).T

    # use nearest neighbor with the tree already created
    tt0 = time()
    dp = DD[xyzT.query(xyzs, n_jobs=-1)[1]] # n_jobs=-1 to use all available cores
    print('NP = %10d: %0.4f sec' % (NP, time() - tt0))
    
if do_plot:
    # plotting
    plt.close('all')
    fig = plt.figure()

    ax = fig.add_subplot(111)
    ax.plot(zp,dp,'.b')
    ax.grid(True)
    ax.axis('square')
    ax.plot([-1,0], [-1,0], '-r', linewidth=2)

    plt.show()

