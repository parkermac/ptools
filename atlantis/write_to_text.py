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
z_dict = {0:0, 1:-5, 2:-25, 3:-50, 4:-100, 5:-150, 6:-350}


#%% setup input locations

in_dir0 = Ldir['parent'] + 'ptools_output/atlantis/'

in_dir = in_dir0 + 'means_and_fluxes/'

dd = os.listdir(in_dir)

# these should be sorted by date
stvwa_list = [item for item in dd if 'stvwa' in item]
trans_list = [item for item in dd if 'trans' in item]

dt0 = datetime(2006,1,1)

fid1 = open(in_dir + 'salish_avg.dat', 'wt')
fid1.write('%5s %5s %5s %15s %15s %15s %15s %15s\n' %
    ('Time', 'nPoly', 'Depth', 'w', 'Temp', 'Salinity', 'Area', 'Volume'))
fid1.write('%5s %5s %5s %15s %15s %15s %15s %15s\n' %
    ('[12h]', '[#]', '[m]', '[1e-6 m/s]', '[degC]', '[psu]', '[km2]', '[km3]'))
for item in stvwa_list:
    dt = datetime.strptime(item[1:11], '%Y.%m.%d')
    nd = (dt - dt0).days
    timestep = nd*2 + 1
    # The file at noon of 2006.01.01 has timestep = 1
    # and the file at noon the next day has timestep = 3
    # because timesteps are in 12 hour increments.
    stvwa_dict = pickle.load(open(in_dir + item, 'rb'))
    for npoly in stvwa_dict.keys():
        stvwa_arr = stvwa_dict[npoly]
        for iz in range(stvwa_arr.shape[0]):
            depth = -z_dict[iz+1]
            salt = stvwa_arr[iz, 0]
            temp = stvwa_arr[iz, 1]
            vol = stvwa_arr[iz, 2]/1e9
            w = stvwa_arr[iz, 3]*1e6
            area = stvwa_arr[iz, 4]/1e6
            fid1.write('%5d %5d %5d %15.7e %15.7e %15.7e %15.7e %15.7e\n' %
                (timestep, npoly, depth, w, temp, salt, area, vol))
fid1.close()

fid2 = open(in_dir + 'salish_flux.dat', 'wt')
fid2.write('%5s %5s %5s %5s %15s\n' %
    ('Time', 'nPoly', 'Face', 'Depth', 'Flux'))
fid2.write('%5s %5s %5s %5s %15s\n' %
    ('[12h]', '[#]', '[#]', '[m]', '[m3/s]'))
for item in trans_list:
    dt = datetime.strptime(item[1:11], '%Y.%m.%d')
    nd = (dt - dt0).days
    timestep = nd*2 + 1
    trans_dict = pickle.load(open(in_dir + item, 'rb'))
    for npoly in trans_dict.keys():
        trans_arr = trans_dict[npoly]
        for iz in range(trans_arr.shape[0]):
            depth = -z_dict[iz+1]
            for face in range(trans_arr.shape[1]):
                flux = trans_arr[iz, face]
                fid2.write('%5d %5d %5d %5d %15.7e\n' %
                    (timestep, npoly, face, depth, flux))
fid2.close()


