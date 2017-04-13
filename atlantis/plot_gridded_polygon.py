# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 08:57:08 2016

@author: PM5

Code to plot our polygons, and determine shared faces.

We should move the shared faces code over to preprocessing
in polygon_on_grid.py.

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

pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import pickle
import numpy as np

save_plots = False
# set this to True to save plot png's to a folder
# set to False to show on screen

#%% setup input location
in_dir0 = Ldir['parent'] + 'ptools_output/atlantis/'
in_dir = in_dir0 + 'gridded_polygons/'

if save_plots:
    out_dir = in_dir + 'plots/'
    Lfun.make_dir(out_dir, clean=True)

#%% get ROMS history file

R_in_dir0 = Ldir['parent'] + 'roms/output/salish_2006_4_lp/'
dt = datetime(2006,7,29)
f_string = 'f' + dt.strftime('%Y.%m.%d')
R_in_dir = R_in_dir0 + f_string + '/'
R_fn = R_in_dir + 'low_passed.nc'
G = zrfun.get_basic_info(R_fn, only_G=True)
lon = G['lon_psi']
lat = G['lat_psi']
# get the vectors that describe the model plaid grid
Lon = lon[0,:].flatten()
Lat = lat[:,0].flatten()

#%% load Atlantis polygon info into a DataFrame

pfn = (Ldir['parent'] + 'PROJECTS/LLTK/Atlantis/Puget_Sound_HydroAtlantis/' +
        'AtlantisBoxInfo_toParker.xlsx')
df = pd.read_excel(pfn, sheetname='BoxVertices')

#%% load polygon results

gpoly_dict = pickle.load(open(in_dir + 'gpoly_dict.p', 'rb'))

#%% plotting

#plt.close('all')

counter = 0
for npoly in [4, 18]:#gpoly_dict.keys():

    this_poly = df[df.box_id==npoly]
    lon_poly = this_poly.Long.values
    lat_poly = this_poly.Lat.values
    
    per_dict = gpoly_dict[npoly]['per_dict']
    ji_rho_in = gpoly_dict[npoly]['ji_rho_in']

    i_dict = gpoly_dict[npoly]['i_dict']
    j_dict = gpoly_dict[npoly]['j_dict']

    if counter == 0:
        fig = plt.figure(figsize=(16, 10))

        # set up panel 1
    
        ax1 = fig.add_subplot(121)
        pfun.add_coast(ax1)
        ax1.axis([-124, -122, 46.8, 49.2])
        pfun.dar(ax1)
        ax1.set_xlabel('Longitude')
        ax1.set_ylabel('Latitude')
    # plot the original polygon, with a star at the start and
    # a segment indicating direction
    ax1.plot(lon_poly, lat_poly,'-r')
    ax1.plot(lon_poly[0], lat_poly[0],'-*m', markersize=12)
    ax1.plot(lon_poly[:2], lat_poly[:2],'-m', linewidth=3)

    #ax1.set_title('Polygon ' + str(npoly))

    # set up panel 2
    if counter == 0:
        ax2 = fig.add_subplot(122)
        pfun.add_coast(ax2)
        ax2.axis([lon_poly.min() - .1, lon_poly.max() + .1,
                 lat_poly.min() - .05, lat_poly.max() + .05,])
        pfun.dar(ax2)
        ax2.set_xlabel('Longitude')
        ax2.set_ylabel('Latitude')
    # plot the original polygon, with a star at the start and
    # a segment indicating direction
    ax2.plot(lon_poly, lat_poly,'-r')
    ax2.plot(lon_poly[0], lat_poly[0],'-*m', markersize=12)
    ax2.plot(lon_poly[:2], lat_poly[:2],'-m', linewidth=3)
    # plot gridded perimeter
    for iseg in i_dict.keys():
        sLon = Lon[i_dict[iseg]]
        sLat = Lat[j_dict[iseg]]
        # plot the points on the psi grid for each segment (stairstep lines)
        ax2.plot(sLon, sLat, '-xb')
        # starting points
        ax2.plot(sLon[0], sLat[0], '-*b', markersize=12)
        # ending points
        ax2.plot(sLon[-1], sLat[-1], '-oy', markersize=12, alpha=.2)
    # plot sign of inflow at u- or v-grid points
    for iseg in per_dict.keys():
        seg = per_dict[iseg]
        if len(seg) > 0:
            for nn in range(seg.shape[0]):
                if seg[nn,2] == 0: # u_grid
                    ax2.text(G['lon_u'][seg[nn,0],seg[nn,1]],
                            G['lat_u'][seg[nn,0],seg[nn,1]],
                            str(seg[nn,3]),
                            horizontalalignment='center',
                            verticalalignment='center')
                elif seg[nn,2] == 1: # v_grid
                    ax2.text(G['lon_v'][seg[nn,0],seg[nn,1]],
                            G['lat_v'][seg[nn,0],seg[nn,1]],
                            str(seg[nn,3]),
                            horizontalalignment='center',
                            verticalalignment='center')
        else:
            pass
    # plot points on the rho grid that are inside the polygon
    ax2.plot(G['lon_rho'][ji_rho_in[:,0], ji_rho_in[:,1]],
              G['lat_rho'][ji_rho_in[:,0], ji_rho_in[:,1]], 'go')
              
    counter += 1

    if save_plots:
        plt.savefig(out_dir + 'poly_' + str(npoly) + '.png')
        plt.close()
    else:
        plt.show()