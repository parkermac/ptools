#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Further process data by region, binned by month and depth.  At this point we merge with the Ecology data.

"""

# imports
import pandas as pd
import numpy as np
import pickle

# directory for processed data to use
dir1 = '../../ptools_output/collias/'
dir2 = '../../ptools_output/ecology/'

# variables to process
col_list = ['Temp. (deg C)', 'Salinity', 'DO (mg L-1)', 'NO3 (uM)']

# years to process
year_list = range(1932, 2018)

A_dict = dict()
for region in range(1,9):
    print('working on region ' + str(region))
    a = pd.read_pickle(dir1 + 'region_' + str(region) + '.p')
    a = a.reindex(columns=(['Z (m)'] + col_list))
    a_ecy = pd.read_pickle(dir2 + 'region_' + str(region) + '.p')
    a_ecy = a_ecy.reindex(columns=(['Z (m)'] + col_list))
    # a has Date as an index
    a = pd.concat((a,a_ecy))

    # get mean values at a specific depth, binned by month
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
            
    A0.index.name = 'Date'
    A10.index.name = 'Date'
    A30.index.name = 'Date'
            
    A_dict[region] = (A0, A10, A30)
    
out_fn = dir1 + 'region_month_z_means.p'
pickle.dump(A_dict, open(out_fn, 'wb'))
