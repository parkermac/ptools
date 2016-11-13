"""
Code to plot all the polygons. Do they really not share faces?

RESULT: Yes they do.

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


#%% setup input location
in_dir0 = Ldir['parent'] + 'ptools_output/atlantis/'
in_dir = in_dir0 + 'gridded_polygons/'

#%% load Atlantis polygon info into a DataFrame

pfn = (Ldir['parent'] + 'PROJECTS/LLTK/Atlantis/Puget_Sound_HydroAtlantis/' +
        'AtlantisBoxInfo_toParker.xlsx')
df = pd.read_excel(pfn, sheetname='BoxVertices')

#%% load polygon results

gpoly_dict = pickle.load(open(in_dir + 'gpoly_dict.p', 'rb'))

#%% plotting

fig = plt.figure(figsize=(16, 10))

# set up panel 1
ax = fig.add_subplot(111)
pfun.add_coast(ax)
ax.axis([-124, -122, 46.8, 49.2])
pfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

for npoly in gpoly_dict.keys():
    
    per_dict = gpoly_dict[npoly]['per_dict']
    NFACE_ALT = len(per_dict)                
    

    this_poly = df[df.box_id==npoly]
    lon_poly = this_poly.Long.values
    lat_poly = this_poly.Lat.values

    # plot the original polygon, with a star at the start and
    # a segment indicating direction
    lh = ax.plot(lon_poly, lat_poly,'*-')
    
    NFACE = len(lon_poly)-1
    for nface in range(NFACE):
        dd = np.random.rand()/300
        x0 = lon_poly[nface]
        x1 = lon_poly[nface+1]
        y0 = lat_poly[nface]
        y1 = lat_poly[nface+1]
        ax.text(np.mean([x0,x1])+dd, np.mean([y0,y1])+dd,
            str(nface), color = lh[0].get_color(), fontsize=10, fontweight='bold')
   
    ax.text(lon_poly.mean(), lat_poly.mean(),
        str(npoly), color = lh[0].get_color(), fontsize=18, fontweight='bold')

    if NFACE != NFACE_ALT:
        print('npoly=%d NFACE=%d NFACE_ALT=%d' % (npoly, NFACE, NFACE_ALT))

plt.show()