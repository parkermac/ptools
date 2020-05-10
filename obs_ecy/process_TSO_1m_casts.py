"""
Process temperature, salinity and DO from Ecology and EC CTD casts.
Uses all depths in 1 m increments.

Good performance - under 30 sec per year.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import netCDF4 as nc4
from time import time

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart()
import zfun

# ***** User Edits

Ldir['gtagex'] = 'cas6_v3_lo8b'

testing = False

# ***** End User Edits

# +++ load Ecology and Canadian CTD cast data +++
dir0 = Ldir['parent'] + 'ptools_data/ecology/'
dir1 = Ldir['parent'] + 'ptools_data/canada/'
sta_df = pd.read_pickle(dir0 + 'sta_df.p')
sta_df_ca = pd.read_pickle(dir1 + 'sta_df.p')
sta_df = pd.concat((sta_df, sta_df_ca), sort=False)

for year in [2017, 2018, 2019]:
    
    print('\n***** ' + str(year) + ' *****')

    # output location
    out_dir = Ldir['parent'] + 'ptools_output/ecology/'
    out_fn = 'TSO_1m_casts_' + Ldir['gtagex'] + '_' + str(year) + '.p'

    #
    Casts = pd.read_pickle(dir0 + 'Casts_' + str(year) + '.p')
    try:
        Casts_ca = pd.read_pickle(dir1 + 'Casts_' + str(year) + '.p')
        Casts = pd.concat((Casts, Casts_ca), sort=False)
    except FileNotFoundError:
        pass

    if testing==True:
        sta_list = [s for s in sta_df.index if 'PSB003' in s]
    else:
        sta_list = [s for s in sta_df.index]

    TSO = pd.DataFrame() # DataFrame to hold all cast data

    tt0 = time()
    # loop over all stations
    for station in sta_list:
        print(station)
        casts = Casts[Casts['Station'] == station]
        casts = casts.set_index('Date')
        # identify a single cast by its DATE
        calldates = casts.index
        dates = calldates.unique() # a short list of unique dates (1 per cast)
    
        # some useful dictionaries for renaming and rescaling
        cast_vn_dict = {'Salinity': 'Salinity', 'Temperature': 'Temp. (deg C)',
            'DO': 'DO (mg L-1)'}
        mod_fac_dict = {'salt': 1, 'temp': 1, 'oxygen': 0.032}
        # model fields have these units AFTER multiplication by mod_fac_dict
        mod_vn_dict = {'salt': 'Mod Salinity', 'temp': 'Mod Temp. (deg C)', 'oxygen': 'Mod DO (mg L-1)'}
    
        # rename the columns
        casts = casts.rename(columns=cast_vn_dict)
    
        # limit the columns
        casts = casts[['Station', 'Z'] + list(cast_vn_dict.values())]
    
        # Set the depths to look at.
        Z_vec = np.arange(0.,-200,-1)

        # loop over all casts at this station
        dates = dates[dates.year==year]
        for dd in dates:

            #Initialize a DataFrame for this station/date
            bc_columns = ['Station', 'Date', 'Z',
                'Salinity', 'Mod Salinity',
                'Temp. (deg C)', 'Mod Temp. (deg C)',
                 'DO (mg L-1)', 'Mod DO (mg L-1)']
            bc = pd.DataFrame(index=Z_vec, columns=bc_columns)
            bc['Date'] = dd
            bc['Station'] = station
            bc['Z'] = Z_vec # we save this as a column and as the index, because
                    # later we will need it after we drop the index during
                    # concatenation into TSO.
    
            # NOTE the brackets around [dd] keep the result as a DataFrame even if
            # we are only pulling out a single row.
            ca = casts.loc[[dd],:]
            ca = ca.set_index('Z')
            for z in ca.index: 
                if z not in Z_vec: 
                    ca.drop(z, inplace=True)
                
            for vn in cast_vn_dict.values():
                cv = ca[vn]
                bc.loc[ca.index,vn] = cv
        
            # finally get the model fields for this day
            date_string = dd.strftime('%Y.%m.%d')
            print('  ' + date_string)
            dir0m = Ldir['LOo'] + 'cast/'
            fnm = dir0m + Ldir['gtagex'] + '/' + station + '_' + date_string + '.nc'
            ds = nc4.Dataset(fnm)
            z = ds['z_rho'][:].squeeze()
            z = z - z[-1] # reference to SSH, so top value is always at 0
        
            i0, i1, fr = zfun.get_interpolant(Z_vec, z, extrap_nan=True)
            for vn in mod_vn_dict.keys():
                mcv = ds[vn][:].squeeze()
                # note we also apply the scaling factor here
                val = ((1-fr)*mcv[i0] + fr*mcv[i1])*mod_fac_dict[vn]
                if isinstance(val, np.ma.MaskedArray):
                    val = val.data
                bc.loc[Z_vec,mod_vn_dict[vn]] = val
            # NOTE: could use bc = bc.dropna() to eliminate all columns that
            # have missing data in them.  I think I will avoid this in case there
            # are casts with good T and s but no DO.  It will not affect the eventual
            # calculation of statistics.
            ds.close()
            TSO = pd.concat((TSO,bc), ignore_index=True, sort=False)
    print('Time to process all stations = %0.1f sec' % (time()-tt0))

    # save output
    TSO.to_pickle(out_dir + out_fn)
