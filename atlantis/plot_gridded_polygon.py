# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 08:57:08 2016

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
import zrfun

pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import pickle

#%% setup input location
in_dir0 = Ldir['parent'] + 'ptools_output/atlantis/'
in_dir = in_dir0 + 'gridded_polygons/'

out_dir = in_dir + 'plots/'
Lfun.make_dir(out_dir, clean=True)


#%% get ROMS history file

R_in_dir0 = Ldir['parent'] + 'roms/output/salish_2006_4_lp/'
dt = datetime(2006,7,29)
f_string = 'f' + dt.strftime('%Y.%m.%d')
R_in_dir = R_in_dir0 + f_string + '/'
R_fn = R_in_dir + 'low_passed.nc'
[G] = zrfun.get_basic_info(R_fn, getS=False, getT=False)
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

plt.close('all')

for npoly in gpoly_dict.keys():

    this_poly = df[df.box_id==npoly]
    lon_poly = this_poly.Long.values
    lat_poly = this_poly.Lat.values

    per_dict = gpoly_dict[npoly]['per_dict']
    ji_rho_in = gpoly_dict[npoly]['ji_rho_in']

    i_dict = gpoly_dict[npoly]['i_dict']
    j_dict = gpoly_dict[npoly]['j_dict']

    fig = plt.figure(figsize=(16, 10))

    # set up panel 1
    ax = fig.add_subplot(121)
    pfun.add_coast(ax)
    ax.axis([-124, -122, 46.8, 49.2])
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # plot the original polygon, with a star at the start and
    # a segment indicating direction
    ax.plot(lon_poly, lat_poly,'-r')
    ax.plot(lon_poly[0], lat_poly[0],'-*m', markersize=12)
    ax.plot(lon_poly[:2], lat_poly[:2],'-m', linewidth=3)

    ax.set_title('Polygon ' + str(npoly))

    # set up panel 2
    ax = fig.add_subplot(122)
    pfun.add_coast(ax)
    ax.axis([lon_poly.min() - .1, lon_poly.max() + .1,
             lat_poly.min() - .05, lat_poly.max() + .05,])
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # plot the original polygon, with a star at the start and
    # a segment indicating direction
    ax.plot(lon_poly, lat_poly,'-r')
    ax.plot(lon_poly[0], lat_poly[0],'-*m', markersize=12)
    ax.plot(lon_poly[:2], lat_poly[:2],'-m', linewidth=3)
    # plot gridded perimeter
    for iseg in i_dict.keys():
        sLon = Lon[i_dict[iseg]]
        sLat = Lat[j_dict[iseg]]
        # plot the points on the psi grid for each segment (stairstep lines)
        ax.plot(sLon, sLat, '-xb')
        # starting points
        ax.plot(sLon[0], sLat[0], '-*b', markersize=12)
        # ending points
        ax.plot(sLon[-1], sLat[-1], '-oy', markersize=12, alpha=.2)
    # plot sign of inflow at u- or v-grid points
    for iseg in per_dict.keys():
        seg = per_dict[iseg]
        if len(seg) > 0:
            for nn in range(seg.shape[0]):
                if seg[nn,2] == 0: # u_grid
                    ax.text(G['lon_u'][seg[nn,0],seg[nn,1]],
                            G['lat_u'][seg[nn,0],seg[nn,1]],
                            str(seg[nn,3]),
                            horizontalalignment='center',
                            verticalalignment='center')
                elif seg[nn,2] == 1: # v_grid
                    ax.text(G['lon_v'][seg[nn,0],seg[nn,1]],
                            G['lat_v'][seg[nn,0],seg[nn,1]],
                            str(seg[nn,3]),
                            horizontalalignment='center',
                            verticalalignment='center')
        else:
            pass
    # plot points on the rho grid that are inside the polygon
    ax.plot(G['lon_rho'][ji_rho_in[:,0], ji_rho_in[:,1]],
              G['lat_rho'][ji_rho_in[:,0], ji_rho_in[:,1]], 'go')


    plt.savefig(out_dir + 'poly_' + str(npoly) + '.png')
    plt.close()