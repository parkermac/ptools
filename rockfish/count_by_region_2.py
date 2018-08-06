#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 14:08:58 2018

@author: pm7

Count the number of particles that end up in a number
of basins, for all the rockfish experiments.

It makes 1829 rows of data:
59 release days [0-58] x 31 settling days [90-120] = 1829
for each experiment.

PERFORMANCE: takes about 16 seconds per experiment, and there are 33
experiments, so this takes a total of 10 minutes.  Not bad.

7/20/2018 Added more detailed information about which polygon
the "lost to land" ones ended up in.

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
outdir = indir + 'counted_by_region_3/'
Lfun.make_dir(outdir, clean=True)
ex_list_raw = os.listdir(indir)
ex_list = []
for ex in ex_list_raw:
    if ('.nc' in ex) and ('grid' not in ex):
        ex_list.append(ex)
ex_list.sort()

if testing:
    ex_list = [ex_list[10]]

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
    
poly_lost_list = []
for po in poly_list:
    poly_lost_list.append('lost_to_land_'+po)
    
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
    # df = pd.DataFrame(columns=['rel_day', 'age'] + poly_list +
    #     ['total_found_in_ploygons','lost_to_ocean','lost_to_land','found_plus_lost'])
    df = pd.DataFrame(columns=['rel_day', 'age'] + poly_list +
        ['total_found_in_ploygons','lost_to_ocean', 'lost_to_land'] + poly_lost_list + ['found_plus_lost'])
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
    # Age gets to 179 for particle 0 and 121 for the last particle.
    # An irregular number of particles are released each day.
    # (see ptools_data/rockfish/rockfish_partruition.csv)
    
    # 1. Find release day of all particles
    nzero = (Age==0).sum(axis=0) - 24
    release_day = nzero/24
    
    rd_unique = np.arange(59)
    # there are only 59 release days, although I expected 60
    Release_day = release_day * np.ones(shape=(NT,1))
    
    # age limits to consider
    age_range = range(90,121)
    
    if testing:
        age_range = [115]
                   
    for age in age_range:

        age_ind = (Age < age).sum(axis=0)
                   
        # and find particle locations at that time
        # using fancy indexing
        lon0 = Lon[age_ind, np.arange(NP)]
        lat0 = Lat[age_ind, np.arange(NP)]
        h0 = H[age_ind, np.arange(NP)]
        rd0 = Release_day[age_ind, np.arange(NP)]

        for rd in rd_unique:
            # loop over all release days
            
            lon1 = lon0[rd0==rd]
            lat1 = lat0[rd0==rd]
            h1 = h0[rd0==rd]
            rd1 = rd0[rd0==rd]
            
            # things for naming
            rd_str = ('0' + str(int(rd)))[-2:]
            age_str = ('00' + str(int(age)))[-3:]
            istr = 'rd' + rd_str + '_age' + age_str
    
            # make masks: True for what we want to keep
            mask_h = h1 >= 50 # at this age they are deeper than 50 m
            # and ones that are shallower are stuck on land
            mask_west = lon1 > -126.5
            mask_south = lat1 > 45.5
            
            mask = mask_h & mask_west & mask_south
            
            lon2 = lon1[mask]
            lat2 = lat1[mask]
            rd2 = rd1[mask]

            lon3 = lon1[~mask_h]
            lat3 = lat1[~mask_h]


            # keep track of how many we lost in the masking
            df.loc[istr,'lost_to_ocean'] = (~mask_west).sum() + (~mask_south).sum()
            df.loc[istr,'lost_to_land'] = (~mask_h).sum()
        
            xy2 = np.stack((lon2,lat2), axis=1)
            xy3 = np.stack((lon3,lat3), axis=1)
        
            # count the number of particles that ended up in each polygon
            # in age categories
            for pn in poly_list:
                p = p_dict[pn]
                p_in = p.contains_points(xy2) # boolean
                df.loc[istr, pn] = p_in.sum()
                
                p_lost = p.contains_points(xy3) # boolean
                df.loc[istr, 'lost_to_land_'+pn] = p_lost.sum()
            
                # and keep track of what we masked away
                # NOTE: if you just look at one age, then the sum of
                # 'ocean' and 'land' (one value each per age)
                # plus the sum of 'total' over all release days
                # then you should recover 10,000, although in practice
                # a few are missing because they are in the gaps between
                # the polygons.  I got 9995 for age = 115 in exp 3_1.
    
    # add columns for release day and age, for convenience.
    for istr in df.index:
        rd = istr[2:4] 
        age = istr[-3:]        
        df.loc[istr, 'rel_day'] = int(rd)
        df.loc[istr, 'age'] = int(age)
            
    df.index.name = 'Filter'
    df['total_found_in_ploygons'] = df[poly_list].sum(axis=1)
    
    # sort in a logical way
    df = df.sort_values(by=['rel_day','age'])
    
    df['found_plus_lost'] = df[['total_found_in_ploygons','lost_to_ocean','lost_to_land']].sum(axis=1)
    
    # and export to csv
    df.to_csv(outdir + 'Exp_' + ex_short + '.csv')
    
    print('  took %0.1f seconds' % (time.time() - tt0))
    


