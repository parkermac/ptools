"""
Process data from WA Dept. of Ecology CTD casts, and
save the results for future use.

"""

import pandas as pd
import numpy as np
import pickle

# where the data is, and where results will be stored
dir0 = '../../ptools_data/ecology/'

year = 2017

if year == 2017:
    ctd_fn = dir0 + 'ParkerMacCready2017CTDDataFeb2018.xlsx'
    sheet_name = '2017Provisional_CTDResults'
    sta_fn = dir0 + 'ParkerMacCreadyCoreStationInfoFeb2018.xlsx'

# PROCESS STATION LOCATION INFO
sta_df = pd.read_excel(sta_fn)
sta_df = sta_df.set_index('Station')
# get locations in decimal degrees
for sta in sta_df.index:
    lat_str = sta_df.loc[sta, 'Lat_NAD83 (deg / dec_min)']
    lat_deg = float(lat_str.split()[0]) + float(lat_str.split()[1])/60
    sta_df.loc[sta,'Latitude'] = lat_deg
    #
    lon_str = sta_df.loc[sta, 'Long_NAD83 (deg / dec_min)']
    lon_deg = float(lon_str.split()[0]) + float(lon_str.split()[1])/60
    sta_df.loc[sta,'Longitude'] = -lon_deg    
sta_df.pop('Lat_NAD83 (deg / dec_min)')
sta_df.pop('Long_NAD83 (deg / dec_min)')
# save result
sta_df.to_pickle(dir0 + 'sta_df_' + str(year) + '.p')

# PROCESS CTD DATA

# data long names; we retain only these fields
data_long_names = ['Salinity', 'Temp', 'Density',
                   'Chla_adjusted', 'DO_raw',
                   'Turbidity', 'Z']

# read in the data (all stations, all casts)
all_casts = pd.read_excel(ctd_fn, sheet_name=sheet_name,
    parse_dates = ['Date'])
    
Casts = pd.DataFrame()

for station in sta_df.index:
    
    print(' - gathering: ' + station)
    casts = all_casts[all_casts['Station'] == station]
    casts = casts.set_index('Date')
    casts['Z'] = -casts['Depth'] # and make a Z column
    casts = casts[data_long_names] # keep only selected columns

    # identify a single cast by its DATE
    alldates = casts.index
    castdates = alldates.unique() # a short list of unique dates (1 per cast)

    # process time series for this station
    for cdate in castdates:
        cast = casts[casts.index==cdate]
        # drop repeat values (aiming for the first of a depth pair)
        zdf = np.diff(cast['Z'])
        zdf = np.concatenate((np.array([1.,]),zdf))
        mask = zdf != 0
        cast = cast[mask]
        cast.sort_values('Z', ascending=False)
        cast = cast[:-5] # drop bottom values (sometimes bad)
        
        cast['Station'] = station
        cast['Date'] = cdate
        
        Casts = pd.concat((Casts, cast), ignore_index=True)
# save result
Casts.to_pickle(dir0 + 'Casts_' + str(year) + '.p')
