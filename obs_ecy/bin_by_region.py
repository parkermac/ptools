#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Processes data from Department of Ecology, combining data from CTD casts
and bottles at the same station.  Keeps data only at 0, 10, and 30 m depth.

Meant to replicate the program by the same name in obs_collias.

"""

# imports
import pandas as pd
import numpy as np

# SSMSP import
import os
import sys
pth = os.path.abspath('../ssmsp')
if pth not in sys.path:
    sys.path.append(pth)
import sfun
import pfun

# +++ load station info +++
dir0 = '../../ptools_data/ecology/'
# load processed station info and data
sta_df = pd.read_pickle(dir0 + 'sta_df.p')

testing = False
if testing:
    year_list = [2006]
    region_list = [3]
else:
    year_list = range(1999, 2018)
    region_list = range(1,9)

# create dicts that associate a region number (Collias system) with a name
# and with a list of modern Ecology station letters
collias_region_names = {1:'Strait of Juan de Fuca',
        2:'Admiralty Inlet',
        3:'Puget Sound Basin',
        4:'Southern Puget Sound',
        5:'Hood Canal',
        6:'Whidbey Basin',
        7:'North Sound',
        8:'San Juan Island Passages'}
ecology_region_lists = {1:['SJF'],
        2:['ADM', 'PTH'],
        3:['PSB', 'ELB', 'EAP', 'CMB'],
        4:['OAK', 'BUD', 'CSE', 'DNA', 'CRR', 'GOR', 'NSQ'],
        5:['HCB'],
        6:['SKG', 'SAR', 'PSS'],
        7:['GRG', 'BLL'],
        8:['RSR']}

for region in region_list:
    print('Working on Region ' + str(region))
    Bc = pd.DataFrame() # DataFrame to hold all bottle and cast data
    
    for year in year_list:
        print(' ' + str(year))
        
        Casts = pd.read_pickle(dir0 + 'Casts_' + str(year) + '.p')
        
        try:
            Bottles = pd.read_pickle(dir0 + 'Bottles_' + str(year) + '.p')
            do_bot = True
        except FileNotFoundError:
            do_bot = False
    

        # loop over all stations
        sta_list = [sta for sta in sta_df.index if sta[:3] in ecology_region_lists[region]]
        for station in sta_list:
    
            casts = Casts[Casts['Station'] == station]
            casts = casts.set_index('Date')    
            # identify a single cast by its DATE
            calldates = casts.index
            castdates = calldates.unique() # a short list of unique dates (1 per cast)
            if do_bot:
                # also for bottles
                bottles = Bottles[Bottles['Station'] == station]   
                bottles = bottles.set_index('Date')    
                # identify a single cast by its DATE
                balldates = bottles.index
                bottledates = balldates.unique() # a short list of unique dates (1 per cast)
                # all valid dates (really bottledates should be a subset of castdates, right?)
                dates = bottledates.union(castdates)
            else:
                dates = castdates
    
            # some useful dictionaries for renaming and rescaling
            cast_vn_dict = {'Salinity': 'Salinity', 'Temperature': 'Temp. (deg C)',
                'DO': 'DO (mg L-1)', 'Z': 'Z (m)'}
            bot_vn_dict = {'NO3(uM)D': 'NO3 (uM)', 'Znom': 'Z (m)'}
    
            # rename the columns
            casts = casts.rename(columns=cast_vn_dict)
            # limit the columns
            casts = casts[['Station'] + list(cast_vn_dict.values())]
            if do_bot:
                bottles = bottles.rename(columns=bot_vn_dict)
                bottles = bottles[['Station'] + list(bot_vn_dict.values())]
    
            # Set the depths to look at.
            Z_list = [0,-10,-30]

            # loop over all casts at this station
            for dd in dates:

                #Initialize a DataFrame for this station/date
                bc_columns = ['Station', 'Date', 'Z (m)',
                    'Salinity', 'Temp. (deg C)', 'DO (mg L-1)','NO3 (uM)']
                bc = pd.DataFrame(index=Z_list, columns=bc_columns)
                bc['Date'] = dd
                bc['Station'] = station
                bc['Z (m)'] = Z_list # we save this as a column and as the index, because
                        # later we will need it after we drop the index during
                        # concatenation into Bc.
        
                # NOTE the brackets around [dd] keep the result as a DataFrame even if
                # we are only pulling out a single row.
                try:
                    ca = casts.loc[[dd],:]
                    ca = ca.set_index('Z (m)')
                    for Zn in Z_list:
                        try:
                            if Zn==0:
                                # OK to extrapolate at surface because shallowest cast data is
                                # always deeper than 0.
                                exn = False
                            else:
                                # but don't extrapolate beyond the deepest cast data when looking
                                # for a given Zn
                                exn = True
                            i0, i1, fr = sfun.get_interpolant(np.array([Zn]), ca.index.values, extrap_nan=exn)
                            for vn in cast_vn_dict.values():
                                cv = ca[vn]
                                bc.loc[Zn,vn] = (1-fr)*cv.iloc[int(i0)] + fr*cv.iloc[int(i1)]
                        except:
                            pass
                except KeyError:
                    pass
        
                if do_bot:
                    try:
                        bo = bottles.loc[[dd],:]
                        bo = bo.set_index('Z (m)')
                        # for bottles we just put the data at Znom, ignoring the fact that
                        # sometimes it is at a different actual depth.  When I looked at the
                        # data it did not seem like a big problem.
                        for Zn in Z_list:
                            for vn in ['NO3 (uM)']:
                                bv = bo.loc[Zn, vn]
                                bc.loc[Zn,vn] = bv
                    except KeyError:
                        #print('vn=%s Zn=%d bv=%0.3f' % (vn, Zn, bv))
                        pass

                Bc = pd.concat((Bc,bc), ignore_index=True, sort=False)
                
    Bc = Bc.set_index('Date')
    Bc.to_pickle('../../ptools_output/ecology/region_' + str(region) + '.p')
