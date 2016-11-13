# -*- coding: utf-8 -*-
"""
Code to plot different estimates of w at the free surface,
and experiment with an algorithm to balance transports.

"""

#%% Imports

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()
import zrfun

import matplotlib.pyplot as plt
plt.close('all')
pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

import numpy as np
from datetime import datetime, timedelta
import pickle
import netCDF4 as nc
import time

import matplotlib as mpl
import matplotlib.cm as cm

do_plot = True

#%% setup input locations
in_dir0 = Ldir['parent'] + 'ptools_output/atlantis/'
in_dir = in_dir0 + 'gridded_polygons/'
R_in_dir0 = Ldir['parent'] + 'roms/output/salish_2006_4_lp/'
# load polygon results
gpoly_dict = pickle.load(open(in_dir + 'gpoly_dict.p', 'rb'))
shared_faces_dict = pickle.load(open(in_dir + 'shared_faces.p', 'rb'))
# specify file to work on
dt0 = datetime(2006,1,1)
ndays = 191
# 209 is 2006.07.29, 200 is 2006.07.20
dt = dt0 + timedelta(days=ndays)
f_string = 'f' + dt.strftime('%Y.%m.%d')
print('Working on day ' + f_string)
R_in_dir = R_in_dir0 + f_string + '/'
R_fn = R_in_dir + 'low_passed.nc'
ds = nc.Dataset(R_fn)

u = ds['u'][:].squeeze()
v = ds['v'][:].squeeze()
w0 = ds['w'][0,-1,:,:].squeeze()

[G, S] = zrfun.get_basic_info(R_fn, getT=False)
zeta = ds['zeta'][0,:,:]
z_rho, z_w =  zrfun.get_z(G['h'], zeta, S)
DA = G['DX'] * G['DY']
DAm = np.ma.masked_where(zeta.mask, DA)
DZ = np.diff(z_w, axis=0)

DYu = G['DY'][:, :-1]
DXv = G['DX'][:-1, :] + np.diff(G['DX'], axis=0)/2

DZu = DZ[:, :, :-1] + np.diff(DZ, axis=2)/2
DZv = DZ[:, :-1, :] + np.diff(DZ, axis=1)/2

DAHu = DZu * DYu
DAHv = DZv * DXv

if do_plot:
    # Setup for PLOTTING
    wmm = .1 # color limits for w in mm/sec
    norm = mpl.colors.Normalize(vmin=-wmm, vmax=wmm)
    cmap = cm.jet
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    aa = [-124, -122, 46.8, 49.2]
    fig = plt.figure(figsize=(18, 10))

# Initialize arrays to store transport data.
face_trans_dict = dict()
poly_conv_dict = dict()

# calculate convergence in each polygon
counter = 0
NPOLY = len(gpoly_dict)
for npoly in range(NPOLY):#gpoly_dict.keys():

    #print('  npoly = ' + str(npoly))

    # we have two objects associated with a given polygon,
    #  * per_dict[iseg] has arrays of boundary information, and
    #  * ji_rho_in is an array of indices of interior points
    per_dict = gpoly_dict[npoly]['per_dict']
    ji_rho_in = gpoly_dict[npoly]['ji_rho_in']
    
    j_in = ji_rho_in[:,0]
    i_in = ji_rho_in[:,1]
    
    # # Find a list of repeated locations in per_dict, because
    # # these represent places where the perimeter folded back
    # # closely on itself and my algorithm gave overlapping points.
    # # This would be OK except that my algorith cannot, in this case
    # # correctly identify which direction is IN or OUT.  Since the
    # # transports should cancel we can just omit these points from the
    # # flux calculation.
    # # Begin by making a list of all locations, then do (*).
    # all_ind = []
    # for nface in per_dict.keys():
    #     per = per_dict[nface] # array for a given face
    #     if len(per) > 0:
    #         for ii in range(per.shape[0]):
    #             item = tuple(per[ii,:])
    #             all_ind.append(item)
    # # all_ind is an array made of ALL the entries in per_dict. 

    # find fluxes through the faces
    NFACE = len(per_dict)                
    face_trans_arr = np.zeros(NFACE)
    face_area_arr = np.zeros(NFACE)

    for nface in range(NFACE):
        per = per_dict[nface]
        # draft_per = per_dict[nface]
        # per = np.array([], dtype=int).reshape((0,4))
        # # trim repeats out of per (*)
        # for ii in range(draft_per.shape[0]):
        #     row = draft_per[ii,:]
        #     if all_ind.count(tuple(row)) == 1:
        #         per = np.concatenate((per,row.reshape(1,4)))
        #     else:
        #         pass
        #         #print('skipping repeat')
        #         #print(row)
                
        JJ = per[:,0]
        II = per[:,1]
        UV = per[:,2]
        PM = per[:,3]
        # vectors of integers
        JJu = JJ[UV==0]
        IIu = II[UV==0]
        PMu = PM[UV==0]
        JJv = JJ[UV==1]
        IIv = II[UV==1]
        PMv = PM[UV==1]
        # shorter vectors of integers, specific to the u- and v-grids
        
        PMu = PMu.reshape((1, len(PMu)))
        PMv = PMv.reshape((1, len(PMv)))

        this_u = u[:, JJu, IIu]
        this_v = v[:, JJv, IIv]

        draft_DAHu = DAHu[:, JJu, IIu]
        draft_DAHv = DAHv[:, JJv, IIv]
        this_DAHu = np.ma.masked_where(this_u.mask, draft_DAHu)
        this_DAHv = np.ma.masked_where(this_v.mask, draft_DAHv)

         # check lengths
        l1u = this_u.compressed().size
        l2u = this_DAHu.compressed().size
        if l1u != l2u:
            print('Warning: U Result vectors are different lengths')

        l1v = this_v.compressed().size
        l2v = this_DAHv.compressed().size
        if l1v != l2v:
            print('Warning: V Result vectors are different lengths')

        # do the integrals
        if l1u>0:
            area_u = this_DAHu.sum()
            trans_u = (this_u * this_DAHu * PMu).sum()
        else:
            area_u = 0.
            trans_u = 0.
        #
        if l1v>0:
            area_v = this_DAHv.sum()
            trans_v = (this_v * this_DAHv * PMv).sum()
        else:
            area_v = 0.
            trans_v = 0.

        #print('iseg = %3d, trans_u = %15d, trans_v = %15d' % (iseg, trans_u, trans_v))
        face_trans = trans_u + trans_v
        face_trans_arr[nface] = face_trans
        face_area = area_u + area_v
        face_area_arr[nface] = face_area
        
        # store results for later
        face_trans_dict[(npoly, nface)] = (face_trans, face_area)
    
    poly_area = DAm[j_in,i_in].sum()
    net_conv = face_trans_arr.sum()       
    if (poly_area > 0):
        poly_mean_w = net_conv/poly_area
    else:
        poly_mean_w = 0.0
        
    # store results for later
    net_face_area = face_area_arr.sum() 
    poly_conv_dict[npoly] = (net_conv, poly_area, net_face_area, NFACE)

    if do_plot:
        if np.abs(poly_mean_w*1000) > wmm/2:
            print('  ** poly_mean_w = %5.3f mm/s' % (poly_mean_w*1000))        
        if counter == 0:
            # set up the plot   
            axp = fig.add_subplot(121)
            axp.axis(aa)
            pfun.dar(axp)
            axp.set_title('Poly Conv w (mm/s)')    
        # plot markers for all points on the rho grid inside the polygon,
        # colored by vertical velocity at the free surface based on the
        # convergence into the polygon    
        axp.plot(G['lon_rho'][j_in,i_in], G['lat_rho'][j_in,i_in], 'o',
            c=m.to_rgba(poly_mean_w*1000), mew=0.0, markersize=2)
        # mew=0.0 sets the markeredgewidth to 0.
        # markersize=1 seems to show nothing
        # c=... sets the markercolor to an rgba color
        # and in some way this is different than markercolor=...        
        axp.text(G['lon_rho'][j_in,i_in].mean(), G['lat_rho'][j_in,i_in].mean(),
            str(npoly), color='k', fontsize=18, fontweight='bold')
    
    counter += 1
#NPOLY = counter
   
# next try balancing transports

# First check to see that all the face transports match
for npoly in range(NPOLY):
    net_conv, poly_area, net_face_area, NFACE = poly_conv_dict[npoly]
    for nface in range(NFACE):
        try:
            face_trans, face_area = face_trans_dict[(npoly, nface)]
            ipoly, iface = shared_faces_dict[(npoly, nface)]
            facing_trans, facing_area = face_trans_dict[(ipoly, iface)]
            
            if np.abs(face_trans + facing_trans) > .001:
                print('(%2d,%3d):(%2d,%3d):' % (npoly, nface, ipoly, iface))
                print('  delta trans = %10.5f' % (face_trans + facing_trans))
                print('    face_trans=%10.2f facing_trans=%10.2f' % (face_trans, facing_trans))
                print('  delta area = %10.2f' % (face_area - facing_area))
        
        except KeyError:
            # presumably this face does not have a match on another polygon
            pass

if do_plot:
    pfun.add_coast(axp)
    plt.show()
                
ds.close()
