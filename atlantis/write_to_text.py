# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 15:03:29 2016

@author: PM5
"""

#%% Imports

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()

import numpy as np
from datetime import datetime, timedelta
import pickle

# specify z levels
z_dict = {0:5, 1:-5, 2:-25, 3:-50, 4:-100, 5:-150, 6:-350}
NLAY = len(z_dict) - 1

#%% setup input locations
in_dir0 = Ldir['parent'] + 'ptools_output/atlantis/'
in_dir_means = in_dir0 + 'means/'
in_dir_fluxes = in_dir0 + 'fluxes/'
dd = os.listdir(in_dir_means)
# these should be sorted by date
stv_list = [item for item in dd if 'stv' in item]

dt0 = datetime(2006,1,1)

# write averages
fid1 = open(in_dir0 + 'salish_avg.dat', 'wt')    
fid1.write('%s\n' % ('      Time Step        Polygon          Depth       Vertical        Average        Average'))
fid1.write('%s\n' % ('         (12)hr         number          Layer       velocity    Temperature       Salinity'))
fid1.write('%s\n\n' % ('                                          [m]    [10^-6*m/s]      [Celsius] [PartsPer1000]'))     
for item in stv_list:
    f_string = item[:11]
    dt = datetime.strptime(f_string[1:], '%Y.%m.%d')
    nd = (dt - dt0).days
    # The file at noon of 2006.01.01 has timestep = 1
    # and the file at noon the next day has timestep = 3
    # because timesteps are in 12 hour increments.    
    stv_dict = pickle.load(open(in_dir_means + item, 'rb'))    
    face_trans_dict = pickle.load(open(in_dir_fluxes + f_string+'_face_trans.p', 'rb'))
    face_ztrans_dict = pickle.load(open(in_dir_fluxes + f_string+'_face_ztrans.p', 'rb'))
    poly_conv_dict = pickle.load(open(in_dir_fluxes + f_string+'_poly_conv.p', 'rb'))
    wz_dict = pickle.load(open(in_dir_fluxes + f_string+'_wz.p', 'rb'))    
    for timestep in [nd*2 + 1, nd*2 + 2]:    
        NPOLY = len(stv_dict.keys())
        for npoly in range(NPOLY):
            stv_arr = stv_dict[npoly]
            for nlay in range(NLAY):
                depth = -z_dict[nlay+1]
                salt = stv_arr[nlay, 0]
                temp = stv_arr[nlay, 1]
                vol = stv_arr[nlay, 2]/1e9
                # get w and horizontal are the BOTTOM of each layer
                if nlay < NLAY -1:
                    w = wz_dict[npoly, nlay+1][0]*1e6 #
                    area = wz_dict[npoly, nlay+1][1]/1e6 # km2
                else:
                    w = 0.
                    area = 0.
                if (not np.isnan(salt)) and (not np.isnan(temp)):
                    fid1.write('%17.7e %17.7e %17.7e %17.7e %17.7e %17.7e\n' %
                        (timestep, npoly, depth, w, temp, salt))
fid1.close()

# write fluxes
fid2 = open(in_dir0 + 'salish_flux.dat', 'wt')
fid2.write('%s\n' % ('        Polygon           Face      Time Step          Depth           Flux'))
fid2.write('%s\n' % ('         number         number         (12)hr          Layer         [m3/s]'))   
for item in stv_list:
    f_string = item[:11]
    dt = datetime.strptime(f_string[1:], '%Y.%m.%d')
    nd = (dt - dt0).days
    # The file at noon of 2006.01.01 has timestep = 1
    # and the file at noon the next day has timestep = 3
    # because timesteps are in 12 hour increments.    
    stv_dict = pickle.load(open(in_dir_means + item, 'rb'))    
    face_trans_dict = pickle.load(open(in_dir_fluxes + f_string+'_face_trans.p', 'rb'))
    face_ztrans_dict = pickle.load(open(in_dir_fluxes + f_string+'_face_ztrans.p', 'rb'))
    poly_conv_dict = pickle.load(open(in_dir_fluxes + f_string+'_poly_conv.p', 'rb'))
    wz_dict = pickle.load(open(in_dir_fluxes + f_string+'_wz.p', 'rb'))
    for timestep in [nd*2 + 1, nd*2 + 2]:        
        NPOLY = len(stv_dict.keys())
        stv_arr = stv_dict[npoly]        
        for npoly in range(NPOLY):
            net_conv, poly_area, poly_zarea, net_face_area, NFACE = poly_conv_dict[npoly]
            for nface in range(NFACE):
                for nlay in range(NLAY):
                    (face_ztrans, face_zarea) = face_ztrans_dict[(npoly, nface, nlay)]                     
                    depth = -z_dict[nlay + 1]
                    flux = face_ztrans                
                    fid2.write('%14d. %14d. %14d. %14d. %14.3f\n' %
                        (npoly, nface + 1, timestep, nlay + 1, flux))
fid2.close()


