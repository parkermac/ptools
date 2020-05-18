#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process data from the Department of Ecology, combining data from CTD casts
and bottles at the same station.  Also including model output.  4 depths.

"""

import pandas as pd
import numpy as np
import netCDF4 as nc4

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart()
import zfun

# ***** User Edits

# specify which model run to use
Ldir['gtagex'] = 'cas6_v3_lo8b'

testing = False

year_list = [2017, 2018, 2019]

if testing:
    year_list = [2017]

# ***** End User Edits

for year in year_list:
    
    print('\n'+str(year)+30*'-=')
    
    # +++ load ecology CTD cast data +++
    dir0 = Ldir['parent'] + 'ptools_data/ecology/'
    # load processed station info and data
    sta_df = pd.read_pickle(dir0 + 'sta_df.p')
    Casts = pd.read_pickle(dir0 + 'Casts_' + str(year) + '.p')
    try:
        # add Canadian data
        dir1 = Ldir['parent'] + 'ptools_data/canada/'
        # load processed station info and data
        sta_df_ca = pd.read_pickle(dir1 + 'sta_df.p')
        sta_df = pd.concat((sta_df, sta_df_ca), sort=False)
        Casts_ca = pd.read_pickle(dir1 + 'Casts_' + str(year) + '.p')
        Casts = pd.concat((Casts, Casts_ca), sort=False)
    except FileNotFoundError:
        pass

    # +++ load ecology bottle data +++
    try:
        Bottles = pd.read_pickle(dir0 + 'Bottles_' + str(year) + '.p')
        do_bottles = True
    except FileNotFoundError:
        do_bottles = False

    # still need to get Canadian bottle data

    if testing==True:
        sta_list = [s for s in sta_df.index if 'SAR003' in s]
    else:
        sta_list = [s for s in sta_df.index]

    # where to save output
    dir11 = Ldir['parent'] + 'ptools_output/ecology/'
    Lfun.make_dir(dir11)

    Bc = pd.DataFrame() # DataFrame to hold all bottle and cast data for this year

    # loop over all stations
    for station in sta_list:
    
        casts = Casts[Casts['Station'] == station]
        casts = casts.set_index('Date')    
        # identify a single cast by its DATE
        calldates = casts.index
        castdates = calldates.unique() # a short list of unique dates (1 per cast)
        # also for bottles
        if do_bottles:
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
        if year == 2019:
            cast_vn_dict = {'Salinity': 'Salinity', 'Temperature': 'Temp. (deg C)',
                'DO': 'DO (mg L-1)'}
        else:
            cast_vn_dict = {'Salinity': 'Salinity', 'Temperature': 'Temp. (deg C)',
                'Chl': 'Chl (mg m-3)', 'DO': 'DO (mg L-1)'}
        bot_vn_dict = {'DIN': 'DIN (uM)'}
        mod_fac_dict = {'salt': 1, 'temp': 1,
            'NO3': 1, 'phytoplankton': 2.5, 'oxygen': 0.032}
        # model fields have these units AFTER multiplication by mod_fac_dict
        mod_vn_dict = {'salt': 'Mod Salinity', 'temp': 'Mod Temp. (deg C)',
            'NO3': 'Mod DIN (uM)', 'phytoplankton': 'Mod Chl (mg m-3)', 'oxygen': 'Mod DO (mg L-1)'}
    
        # rename the columns
        casts = casts.rename(columns=cast_vn_dict)
        if do_bottles:
            bottles = bottles.rename(columns=bot_vn_dict)
    
        # limit the columns
        casts = casts[['Station', 'Z'] + list(cast_vn_dict.values())]
        if do_bottles:
            bottles = bottles[['Station', 'Z', 'Znom'] + list(bot_vn_dict.values())]
    
        # Set the depths to look at.
        # [0,-10,-30] are the bottle nominal depths, but you can get deeper cast data
        # by adding to this list, e.g. -100
        Z_list = [0,-10,-30, -80]

        # loop over all casts at this station
        dates = dates[dates.year==year]
    
        for dd in dates:

            #Initialize a DataFrame for this station/date
            bc_columns = ['Station', 'Date', 'Znom', 'Z',
                'Salinity', 'Temp. (deg C)','Chl (mg m-3)', 'DO (mg L-1)','DIN (uM)',
                'Mod Salinity', 'Mod Temp. (deg C)', 'Mod Chl (mg m-3)',
                'Mod DO (mg L-1)', 'Mod DIN (uM)']
            bc = pd.DataFrame(index=Z_list, columns=bc_columns)
            bc['Date'] = dd
            bc['Station'] = station
            bc['Znom'] = Z_list # we save this as a column and as the index, because
                    # later we will need it after we drop the index during
                    # concatenation into Bc.
    
            # NOTE the brackets around [dd] keep the result as a DataFrame even if
            # we are only pulling out a single row.
            ca = casts.loc[[dd],:]
            ca = ca.set_index('Z')
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
                    i0, i1, fr = zfun.get_interpolant(np.array([Zn]), ca.index.values, extrap_nan=exn)
                    for vn in cast_vn_dict.values():
                        cv = ca[vn]
                        bc.loc[Zn,vn] = (1-fr)*cv.iloc[int(i0)] + fr*cv.iloc[int(i1)]
                except:
                    pass
    
            if do_bottles:
                try:
                    bo = bottles.loc[[dd],:]
                    bo = bo.set_index('Znom')
                    # for bottles we just put the data at Znom, ignoring the fact that
                    # sometimes it is at a different actual depth.  When I looked at the
                    # data it did not seem like a big problem.
                    for Zn in Z_list:
                        for vn in bot_vn_dict.values():
                            bv = bo.loc[Zn, vn]
                            bc.loc[Zn,vn] = bv
                except KeyError:
                    pass

            # finally get the model fields for this day
            date_string = dd.strftime('%Y.%m.%d')
            dir0m = Ldir['LOo'] + 'cast/'
            fnm = dir0m + Ldir['gtagex'] + '/' + station + '_' + date_string + '.nc'
            if testing == True:
                print(fnm)
            try:
                ds = nc4.Dataset(fnm)
                z = ds['z_rho'][:].squeeze()
                z = z - z[-1] # reference to SSH, so top value is always at 0
                for Zn in Z_list:
                    i0, i1, fr = zfun.get_interpolant(np.array([Zn]), z, extrap_nan=True)
                    for vn in mod_vn_dict.keys():
                        try:
                            mcv = ds[vn][:].squeeze()
                            # note we also apply the scaling factor here
                            val = ((1-fr)*mcv[int(i0)] + fr*mcv[int(i1)])*mod_fac_dict[vn]
                            bc.loc[Zn,mod_vn_dict[vn]] = val
                        except IndexError:
                            pass
                ds.close()
            except OSError:
                if testing == True:
                    print('error reading fnm')
                pass
            Bc = pd.concat((Bc,bc), ignore_index=True, sort=False)
            
        
    # save the output to disk
    if testing == False:
        out_fn = 'ObsMod_' + Ldir['gtagex'] + '_'+ str(year) + '.p'
        print('Saving ' + out_fn)
        Bc.to_pickle(dir11 + out_fn)

