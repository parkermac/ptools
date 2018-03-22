#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 14:08:58 2018

@author: pm7

Count the number of particles that end up in a number
of basins, for all the rockfish experiments.

It makes 1829 rows of data:
59 release days [1-59] x 31 settling days [90-120] = 1829
for each experiment.

PERFORMANCE: takes about 14 seconds per experiment, and there are 33
experiments, so this takes a total of 8 minutes.  Not bad.

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
import time

testing = False

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

# create the list of run files
indir = Ldir['parent'] + 'ptools_output/rockfish/'
outdir = indir + 'counted_by_region_2/'
Lfun.make_dir(outdir, clean=True)
ex_list_raw = os.listdir(indir)
ex_list = []
for ex in ex_list_raw:
    if ('.nc' in ex) and ('grid' not in ex):
        ex_list.append(ex)
ex_list.sort()

if testing:
    ex_list = ex_list[:3] #[ex_list[0]]

# make shorter names for experiments
ex_dict = dict()
for ex in ex_list:
    ex_short = ex.replace('rockfish_2006_Experiment_', '')
    ex_short = ex_short.replace('.nc','')
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
    
for ex in ex_list:
    
    tt0 = time.time()
    
    # set up a DataFrame to store results
    df = pd.DataFrame(columns=['rel_day', 'age'] + poly_list + ['total'])
    ex_short = ex_dict[ex]
    
    print('Working on ' + ex_short)
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

    NT, NP = Lon.shape
    # NT = 4320 = 180*24 or hourly for 180 days [0:180]
    # Age is always the integer day, so it repeats 24 times each day.
    # Age gets to 180 for particle 0 and 121 for the last particle.
    # I think 167 particles are released each day for the first 60 days,
    # so I want to isolate these. Actually there are an irregular number
    # of particles released, between 12 and 575 for the first experiment.
    # Odd.
    
    # 1. Find release day of all particles
    nzero = (Age==0).sum(axis=0)
    release_day = nzero/24
    rd_unique = np.unique(release_day)
    Release_day = release_day * np.ones(shape=(NT,1))
    
    age0 = 90
    age1 = 121
    
    if testing:
        rd_unique = rd_unique[:3]
        age0 = 115
        age1 = 121
                   
    for age in range(age0, age1):

        age_ind = (Age <= age).sum(axis=0) 
                   
        # and find particle locations at that time
        # using fancy indexing
        lon1 = Lon[age_ind, np.arange(NP)]
        lat1 = Lat[age_ind, np.arange(NP)]
        h1 = H[age_ind, np.arange(NP)]
        rd1 = Release_day[age_ind, np.arange(NP)]
    
        # Drop particles that end up deeper than 50 m
        # because the shallow ones are typically stuck on land.
        # Also drop particles that escape the domain.
        mask1 = h1 > 50
        mask2 = lon1 > -126.5
        mask3 = lat1 > 45.5
        mask = mask1 & mask2 & mask3
        lon1 = lon1[mask]
        lat1 = lat1[mask]
        rd1 = rd1[mask]
     
        xy1 = np.stack((lon1,lat1), axis=1)
        
        # count the number of particles that ended up in each polygon
        # in age categories
        for pn in poly_list:
            p = p_dict[pn]
            p_in = p.contains_points(xy1) # boolean
            for rd in rd_unique:
                rd_str = ('0' + str(int(rd)))[-2:]
                age_str = ('00' + str(int(age)))[-3:]
                istr = 'rd' + rd_str + '_age' + age_str
                df.loc[istr, pn] = p_in[rd1==rd].sum()
    
    # add columns for release day and age, for convenience.
    for istr in df.index:
        rd = istr[2:4] 
        age = istr[-3:]        
        df.loc[istr, 'rel_day'] = int(rd)
        df.loc[istr, 'age'] = int(age)
            
    df.index.name = 'Filter'
    df['total'] = df[poly_list].sum(axis=1)
    
    # sort in a logical way
    df = df.sort_values(by=['rel_day','age'])
    
    # and export to csv
    df.to_csv(outdir + 'Exp_' + ex_short + '.csv')
    
    print('  took %0.1f seconds' % (time.time() - tt0))
    


