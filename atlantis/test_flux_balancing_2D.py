# -*- coding: utf-8 -*-
"""
Code to calculate convergence in all polygons
and experiment with an algorithm to balance transports
through the faces.

This code is 2D, meaning we don't consider depth dependence and
instead sum up evelything coming into a polygon.

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
import pandas as pd

import matplotlib as mpl
import matplotlib.cm as cm

do_plot = True

#%% load Atlantis polygon info into a DataFrame
pfn = (Ldir['parent'] + 'PROJECTS/LLTK/Atlantis/Puget_Sound_HydroAtlantis/' +
        'AtlantisBoxInfo_toParker.xlsx')
df = pd.read_excel(pfn, sheetname='BoxVertices')
bi = pd.read_excel(pfn, sheetname='BoxInfo')

#%% setup input locations
in_dir0 = Ldir['parent'] + 'ptools_output/atlantis/'
in_dir = in_dir0 + 'gridded_polygons/'
R_in_dir0 = Ldir['parent'] + 'roms/output/salish_2006_4_lp/'

# load polygon results
gpoly_dict = pickle.load(open(in_dir + 'gpoly_dict.p', 'rb'))
shared_faces_dict = pickle.load(open(in_dir + 'shared_faces.p', 'rb'))

# specify ROMS file to work on
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
# DZ on u and v grids
DZu = DZ[:, :, :-1] + np.diff(DZ, axis=2)/2
DZv = DZ[:, :-1, :] + np.diff(DZ, axis=1)/2
# DX and DY on u and v grids
DYu = G['DY'][:, :-1]
DXv = G['DX'][:-1, :] + np.diff(G['DX'], axis=0)/2
# cell areas for u and v grid box faces
DAHu = DZu * DYu
DAHv = DZv * DXv

if do_plot:
    # Setup for PLOTTING
    wmm = .025 # color limits for w [mm/sec]
    # the next three lines set up to plot points with
    # colors related to w
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
    
    this_poly = df[df.box_id==npoly]
    lon_poly = this_poly.Long.values
    lat_poly = this_poly.Lat.values
    
    this_bi = bi[bi.box_id == npoly]
    boundary = this_bi.Boundary_boxes.values
    island = this_bi.Island_boxes.values

    #print('  npoly = ' + str(npoly))

    # we have two objects associated with a given polygon,
    #  * per_dict[iseg] has arrays of boundary information, and
    #  * ji_rho_in is an array of indices of interior points
    per_dict = gpoly_dict[npoly]['per_dict']
    ji_rho_in = gpoly_dict[npoly]['ji_rho_in']
    
    j_in = ji_rho_in[:,0]
    i_in = ji_rho_in[:,1]
    
    # find fluxes through the faces
    NFACE = len(per_dict)                
    face_trans_arr = np.zeros(NFACE)
    face_area_arr = np.zeros(NFACE)

    for nface in range(NFACE):
        per = per_dict[nface]                
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
            pfun.add_coast(axp)
            pfun.dar(axp)
            axp.set_title('Poly Conv w (mm/s)')    
        # plot filled polygons,
        # colored by vertical velocity at the free surface based on the
        # convergence into the polygon    
        # axp.plot(G['lon_rho'][j_in,i_in], G['lat_rho'][j_in,i_in], 'o',
        #     c=m.to_rgba(poly_mean_w*1000), mew=0.0, markersize=2)
        if island == 0:
            axp.fill(lon_poly, lat_poly,
                c=m.to_rgba(poly_mean_w*1000))
        # mew=0.0 sets the markeredgewidth to 0.
        # markersize=1 seems to show nothing
        # c=... sets the markercolor to an rgba color
        # and in some way this is different than markercolor=...        
        # axp.text(G['lon_rho'][j_in,i_in].mean(), G['lat_rho'][j_in,i_in].mean(),
        #     str(npoly), color='k', fontsize=18, fontweight='bold')
    
    counter += 1

ds.close()

if do_plot:
    plt.show()
                  
# Next try balancing transports

# First check to see that all the face transports match
# RESULT: it works perfectly - all face fluxes match their counterparts
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

# save originals
orig_poly_conv_dict = poly_conv_dict.copy()
orig_face_trans_dict = face_trans_dict.copy()

# Next try to adjust all polygons to have conv = 0
NITER = 100
for iii in range(NITER):
    
    new_poly_conv_dict = poly_conv_dict.copy()
    for npoly in range(NPOLY):
        net_conv, poly_area, net_face_area, NFACE = new_poly_conv_dict[npoly]
        new_poly_conv_dict[npoly] = (0.0, poly_area, net_face_area, NFACE)
    
    new_face_trans_dict = face_trans_dict.copy()

    for npoly in range(NPOLY):
        net_conv, poly_area, net_face_area, NFACE = poly_conv_dict[npoly]
        new_net_conv, poly_area, net_face_area, NFACE = new_poly_conv_dict[npoly]
        dconv = new_net_conv - net_conv
    
        if net_face_area != 0.0:
            dconv_a = dconv / net_face_area
        else:
            dconv_a = 0.0
    
        for nface in range(NFACE):
            face_trans, face_area = face_trans_dict[(npoly, nface)]
            if face_trans != 0.0 and face_area != 0.0:
                new_face_trans_dict[(npoly, nface)] = (face_trans + dconv_a*face_area, face_area)
            else:
                pass # keep original values

    new_face_trans_dict_copy = new_face_trans_dict.copy()
    for npoly in range(NPOLY):
        net_conv, poly_area, net_face_area, NFACE = poly_conv_dict[npoly]
        for nface in range(NFACE):
            try:
                new_face_trans, face_area = new_face_trans_dict[(npoly, nface)]
                ipoly, iface = shared_faces_dict[(npoly, nface)]
                new_facing_trans, facing_area = new_face_trans_dict_copy[(ipoly, iface)]
                fact = (new_face_trans + new_facing_trans)/2
                new_face_trans_dict[(npoly, nface)] = (new_face_trans - fact, face_area)
            except KeyError:
                # presumably this face does not have a match on another polygon
                pass

    counter = 0            
    for npoly in range(NPOLY):
        net_conv, poly_area, net_face_area, NFACE = poly_conv_dict[npoly]
        shelf = []
        for nface in range(NFACE):
            new_face_trans, face_area = new_face_trans_dict[(npoly, nface)]
            shelf.append(new_face_trans)
        new_conv = np.array(shelf).sum()
        new_poly_conv_dict[npoly] = (new_conv, poly_area, net_face_area, NFACE)
    
        
        if do_plot and iii == NITER - 1:
            
            this_poly = df[df.box_id==npoly]
            lon_poly = this_poly.Long.values
            lat_poly = this_poly.Lat.values
            
            this_bi = bi[bi.box_id == npoly]
            boundary = this_bi.Boundary_boxes.values
            island = this_bi.Island_boxes.values
            
        
            if (poly_area > 0):
                poly_mean_w = new_conv/poly_area
            else:
                poly_mean_w = 0.0
                
            if np.abs(poly_mean_w*1000) > wmm/2:
                print('  ** poly_mean_w = %5.3f mm/s' % (poly_mean_w*1000))        
            if counter == 0:
                # set up the plot   
                axp = fig.add_subplot(122)
                axp.axis(aa)
                pfun.add_coast(axp)
                pfun.dar(axp)
                axp.set_title('After iteration')
            if island == 0:   
                axp.fill(lon_poly, lat_poly,
                    c=m.to_rgba(poly_mean_w*1000))
            
            counter += 1
            
    face_trans_dict = new_face_trans_dict.copy()
    poly_conv_dict = new_poly_conv_dict.copy()
    
# look at detailed results
for npoly in range(NPOLY):
    net_conv, poly_area, net_face_area, NFACE = poly_conv_dict[npoly]
    for nface in range(NFACE):
        face_trans, face_area = face_trans_dict[(npoly, nface)]
        orig_face_trans, face_area = orig_face_trans_dict[(npoly, nface)]
        print('(%2d,%3d):(%2d,%3d):' % (npoly, nface, ipoly, iface))
        print('  delta trans = %10.5f' % (face_trans - orig_face_trans))
    
    