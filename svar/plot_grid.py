#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 17 12:57:19 2018

@author: pm7

Code to make a nice figure of the aestus grid.
"""

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart('aestus1', 'A1')
Ldir['gtagex'] = Ldir['gtag'] + '_' + 'ae1'
import zrfun

pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun



dir0 = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[0]=='f']
fn = dir0 + f_list[0] + '/ocean_his_0001.nc'

G = zrfun.get_basic_info(fn, only_G=True)

x = G['lon_psi']
y = G['lat_psi']
h = G['h'][1:-1, 1:-1]
m = G['mask_rho'][1:-1, 1:-1]
h[m==0] = np.nan
z = -h

three = False

xx = G['lon_rho'][1:-1, 1:-1]
yy = G['lat_rho'][1:-1, 1:-1]
hh = G['h'][1:-1, 1:-1]
zz = -hh

# plotting

plt.close('all')

fig = plt.figure(figsize=(11,6))

if three:
    ax = fig.add_subplot(111, projection='3d')
else:
    ax = fig.add_subplot(111)

if three:
    
    ax.plot_surface(xx, yy,zz, rcount=1000, ccount=1000)

else:
    
    cs = ax.pcolormesh(x, y, z, cmap='viridis', vmin=-75, vmax=0)
    pfun.dar(ax)  
    fig.colorbar(cs)
    
    ax.contour(xx,yy,zz, [-70, -60, -50, -40, -30, -20, -10],
               colors='k', linewidths=0.5, linestyles='solid')
    
    ax.set_xlabel('Longitude (deg)')
    ax.set_ylabel('Latitude (deg)')
    ax.set_title('Model Domain and Bathymetry (m)')
    
    ax.set_xlim(-1, 3)
    ax.set_ylim(44, 46)
    
    clr = 'brown'
    ax.plot([0.02, 1.5, 1.5, 0.02, 0.02], [44.9, 44.9, 45.1, 45.1, 44.9],
            ':', color=clr, linewidth=3)
    ax.text(0.2, 45.15, 'Analysis Region', fontsize=14, fontweight='bold',
            color=clr)


plt.show()
