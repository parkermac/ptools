#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 12:11:45 2018

Plot polygons.

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()

pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun
import matplotlib.pyplot as plt
import pickle
import numpy as np

# polygons
poly_dir = Ldir['parent'] + 'ptools_data/rockfish/polygons_rockfish/'
poly_list = ['admiralty', 'hood_canal', 'jdf', 'main_basin', 'outer_coast',
    'san_juans', 'sog', 'south_sound', 'whidbey']
v_dict = dict() # V is for the Vertices of the polygons
for pn in poly_list:
    poly_fn = poly_dir + pn + '.p'
    poly = pickle.load(open(poly_fn, 'rb'))
    px = poly['lon_poly']
    py = poly['lat_poly']
    v = np.ones((len(px),2))
    v[:,0] = px
    v[:,1] = py
    v_dict[pn] = v

# plotting
plt.close('all')

fig = plt.figure(figsize=(12,8))

ax = fig.add_subplot(121)
pfun.add_coast(ax)
ax.set_xlim(-127, -121.5)
ax.set_ylim(43, 50.5)
pfun.dar(ax)
for name in poly_list:
    lon = v_dict[name][:,0]
    lat = v_dict[name][:,1]
    lon = np.append(lon, lon[0])
    lat = np.append(lat, lat[0])
    h = ax.plot(lon, lat, '-', linewidth=2)
    ax.text(lon.mean(), lat.mean(), name, fontsize=10,
            horizontalalignment='center', color=h[0].get_color(),
            fontweight='bold')
 
ax = fig.add_subplot(122)
pfun.add_coast(ax)
ax.set_xlim(-125.5, -122)
ax.set_ylim(47, 49.5)
pfun.dar(ax)
for name in poly_list:
    lon = v_dict[name][:,0]
    lat = v_dict[name][:,1]
    lon = np.append(lon, lon[0])
    lat = np.append(lat, lat[0])
    h = ax.plot(lon, lat, '-', linewidth=2)
    ax.text(lon.mean(), lat.mean(), name, fontsize=10,
            horizontalalignment='center', color=h[0].get_color(),
            fontweight='bold')

plt.show()
    
