#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Look at data by region, binned by month and depth

Regions (the first digit is stored ar Num0 in Bottles):
        100        Strait of Juan de Fuca
        200        Admiralty Inlet
        300        Puget Sound Basin
        400        Southern Puget Sound
        500        Hood Canal
        600        Whidbey Basin
        700        North Sound
        800        San Juan Island Passages

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()

import pandas as pd
import numpy as np
import pickle

dir0 = Ldir['parent'] + 'ptools_data/collias/'
# load processed station info and data
sta_df = pd.read_pickle(dir0 + 'sta_df.p')

# directory for processed data to use
dir1 = Ldir['parent'] + 'ptools_output/collias/'

# variables to process
col_list = ['Temp. (deg C)', 'Salinity',
            'DO (mg L-1)', 'NO3 (mg L-1)', 'NO2 (mg L-1)', 'SiOH4 (mg L-1)']

# years to process
year_list = range(1932, 1976)

A_dict = dict()
for region in range(1,9):
    print('working on region ' + str(region))
    a = pd.read_pickle(dir1 + 'region_' + str(region) + '.p')
    # a has Date as an index

    # get mean values at a specific depth, bined by month
    az = a[a['Z (m)']== 0]
    A0 = pd.DataFrame(columns=col_list)
    for y in year_list:
        for m in range(1,13):
            aa = az[(az.index.year==y) & (az.index.month==m)]
            aa = aa.reindex(columns=col_list) # drops 'Z (m)' column 
            aam = aa.mean(axis=0) # nanmean of all data columns
            # at this month and z level
            A0.loc[pd.datetime(y,m,15),:] = aam
        
    az = a[a['Z (m)']== -10]
    A10 = pd.DataFrame(columns=col_list)
    for y in year_list:
        for m in range(1,13):
            aa = az[(az.index.year==y) & (az.index.month==m)]
            aa = aa.reindex(columns=col_list)
            aam = aa.mean(axis=0)
            A10.loc[pd.datetime(y,m,15),:] = aam
        
    az = a[a['Z (m)']== -30]
    A30 = pd.DataFrame(columns=col_list)
    for y in year_list:
        for m in range(1,13):
            aa = az[(az.index.year==y) & (az.index.month==m)]
            aa = aa.reindex(columns=col_list)
            aam = aa.mean(axis=0)
            A30.loc[pd.datetime(y,m,15),:] = aam
            
    A_dict[region] = (A0, A10, A30)
    
out_fn = dir1 + 'region_month_z_means.p'
pickle.dump(A_dict, open(out_fn, 'wb'))
