# -*- coding: utf-8 -*-
"""
Count the number of particles that end up in a number
of basins, for all the rockfich experiments.
"""

# setup
import os
import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

import numpy as np
import netCDF4 as nc
import matplotlib.path as mpath
import pickle
import pandas as pd

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

# create the list of run files
indir = odir00 = Ldir['parent'] + 'ptools_output/rockfish/'
ex_list_raw = os.listdir(indir)
ex_list = []
for ex in ex_list_raw:
    if ('.nc' in ex) and ('grid' not in ex):
        ex_list.append(ex)
ex_list.sort()

# testing
#ex_list = ex_list[:3]

# make shorter names for experiments
ex_dict = dict()
for ex in ex_list:
    ex_short = ex.replace('rockfish_2006_Experiment', 'ex')
    ex_dict[ex] = ex_short

# polygons
poly_dir = Ldir['parent'] + 'ptools_data/rockfish/polygons_rockfish/'
poly_list = ['admiralty', 'hood_canal', 'jdf', 'main_basin', 'outer_coast',
    'san_juans', 'sog', 'south_sound', 'whidbey']
v_dict = dict() # V is for the Vertices of the polygons
p_dict = dict()
for pn in poly_list:
    poly_fn = poly_dir + pn + '.p'
    poly = pickle.load(open(poly_fn, 'rb'))
    px = poly['lon_poly']
    py = poly['lat_poly']
    v = np.ones((len(px),2))
    v[:,0] = px
    v[:,1] = py
    v_dict[pn] = v
    p_dict[pn] = mpath.Path(v)
    
# set up a DataFrame to store results
df = pd.DataFrame(columns=poly_list)

for ex in ex_list:
    # get data
    ds = nc.Dataset(indir + ex)
    # tracks are stored (time, particle)
    Lon = ds['lon'][:]
    Lat = ds['lat'][:]
    Z = ds['z'][:]
    H = ds['h'][:]
    Age = ds['age'][:]
    ot = ds['ot'][:]
    ds.close()
    # subsample the number of particles
    pstep = 1
    lon = Lon[:,::pstep]
    lat = Lat[:,::pstep]
    z = Z[:,::pstep]
    h = H[:,::pstep]
    age = Age[:,::pstep]
    NT, NP = lon.shape

    # select particles that are 120 days old
    age_ind = (age <= 120).sum(axis=0)
    # and find particle locations at that time
    # using fancy indexing
    lon1 = lon[age_ind, np.arange(NP)]
    lat1 = lat[age_ind, np.arange(NP)]
    xy1 = np.stack((lon1,lat1), axis=1)
    # count the number of particles that ended up in each polygon
    for pn in poly_list:
        p = p_dict[pn]
        p_in = p.contains_points(xy1) # boolean
        #print("%s has %d points" % (pn, p_in.sum()))
        df.loc[ex_dict[ex], pn] = p_in.sum()
    


