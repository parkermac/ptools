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

#%% setup input locations

whichyear = 2005
if Ldir['env'] == 'pm_mac': # mac version
    if whichyear == 2006:
        in_dir0 = Ldir['parent'] + 'ptools_output/atlantis_mac_2006/'
        out_dir0 = Ldir['parent'] + 'ptools_output/atlantis_textfiles_mac_2006/'
    elif whichyear == 2005:
        in_dir0 = Ldir['parent'] + 'ptools_output/atlantis_mac_2005/'
        out_dir0 = Ldir['parent'] + 'ptools_output/atlantis_textfiles_mac_2005/'
elif Ldir['env'] == 'fjord': # fjord version
    if whichyear == 2006:
        in_dir0 = Ldir['parent'] + 'ptools_output/atlantis_fjord_2006/'
        out_dir0 = Ldir['parent'] + 'ptools_output/atlantis_textfiles_fjord_2006/'
    elif whichyear == 2005:
        in_dir0 = Ldir['parent'] + 'ptools_output/atlantis_fjord_2005/'
        out_dir0 = Ldir['parent'] + 'ptools_output/atlantis_textfiles_fjord_2005/'

Lfun.make_dir(out_dir0, clean=True)

# specify z levels
z_dict = {0:5, 1:-5, 2:-25, 3:-50, 4:-100, 5:-150, 6:-350}
NLAY = len(z_dict) - 1

# starting day
dt0 = datetime(whichyear,1,1)

in_dir_means = in_dir0 + 'means/'
in_dir_fluxes = in_dir0 + 'fluxes/'

dd = os.listdir(in_dir_means)
dd.sort()
# these should be sorted by date
stv_list = [item for item in dd if 'stv' in item]

mo_list = range(1,13)

for mo in mo_list:
    
    print(mo)

    this_stv_list = []
    for item in stv_list:
        f_string = item[:11]
        dt = datetime.strptime(f_string[1:], '%Y.%m.%d')
        if dt.month == mo:
            this_stv_list.append(item)
            
    if len(this_stv_list) > 0:
        yr_string = str(dt.year)
        mo_string = ('0' + str(mo))[-2:]
        out_dir = in_dir0 + ('%s_%s/' % (yr_string, mo_string))
        Lfun.make_dir(out_dir, clean=True)
        
        print('writing to ' + out_dir)
        sys.stdout.flush()

        # prepare to organize output
        avg_dict = dict()
        flux_dict = dict()
        timestep_list = []
        nface_dict = dict()
    
        # step through the days in this month
        for item in this_stv_list:
            f_string = item[:11]
            dt = datetime.strptime(f_string[1:], '%Y.%m.%d')
        
            nd = (dt - dt0).days
            # The file at noon of 2006.01.01 has timestep = 1
            # and the file at noon the next day has timestep = 3
            # because timesteps are in 12 hour increments.

            stv_dict = pickle.load(open(in_dir_means + item, 'rb'))    
            face_ztrans_dict = pickle.load(open(in_dir_fluxes + f_string+'_face_ztrans.p', 'rb'))
            poly_conv_dict = pickle.load(open(in_dir_fluxes + f_string+'_poly_conv.p', 'rb'))
            wz_dict = pickle.load(open(in_dir_fluxes + f_string+'_wz.p', 'rb'))

            NPOLY = len(stv_dict.keys())

            for timestep in [nd*2 + 1, nd*2 + 2]:

                timestep_list.append(timestep)
                      
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
                        avg_dict[(timestep, npoly, nlay)] = (depth, w, temp, salt)
                        # salt and temp may still have nan values the will be omitted later
            
                for npoly in range(NPOLY):
                    net_conv, poly_area, poly_zarea, net_face_area, NFACE = poly_conv_dict[npoly]
                    nface_dict[npoly] = NFACE
                    for nface in range(NFACE):
                        for nlay in range(NLAY):                    
                            (face_ztrans, face_zarea) = face_ztrans_dict[(npoly, nface, nlay)]                     
                            flux = face_ztrans                    
                            flux_dict[(npoly, nface, timestep, nlay)] = flux
                              
        # write averages
        fid1 = open(out_dir + 'salish_avg.dat', 'wt')    
        fid1.write('%s\n' % ('      Time Step        Polygon          Depth       Vertical        Average        Average'))
        fid1.write('%s\n' % ('         (12)hr         number          Layer       velocity    Temperature       Salinity'))
        fid1.write('%s\n\n' % ('                                          [m]    [10^-6*m/s]      [Celsius] [PartsPer1000]'))
        for timestep in timestep_list:
            for npoly in range(NPOLY):
                for nlay in range(NLAY):
                    depth, w, temp, salt = avg_dict[(timestep, npoly, nlay)]
                    if (not np.isnan(salt)) and (not np.isnan(temp)):
                        fid1.write('%17.7e %17.7e %17.7e %17.7e %17.7e %17.7e\n' %
                            (timestep, npoly, depth, w, temp, salt))
        fid1.close()

        # write fluxes
        fid2 = open(out_dir + 'salish_flux.dat', 'wt')
        fid2.write('%s\n' % ('        Polygon           Face      Time Step          Depth           Flux'))
        fid2.write('%s\n' % ('         number         number         (12)hr          Layer         [m3/s]'))   
        for npoly in range(NPOLY):
            for nface in range(nface_dict[npoly]):
                for timestep in timestep_list:
                    for nlay in range(NLAY):
                        flux = flux_dict[(npoly, nface, timestep, nlay)]
                        fid2.write('%14d. %14d. %14d. %14d. %14.3f\n' %
                            (npoly, nface+1, timestep, nlay+1, flux))
        fid2.close()



