"""
Process data from WA Dept. of Ecology CTD casts, and
save the results for future use.

"""

import pandas as pd
import numpy as np
import pickle

# where the data is, and where results will be stored
dir0 = '../../ptools_data/ecology/'
sta_df = pd.read_pickle(dir0 + 'sta_df.p')

# PROCESS CTD DATA
read_new = True
for year in [2018]:#range(1999,2018):
    print('\nYEAR = ' + str(year))
    
    
    if year == 2017:
        read_new = True
        ctd_fn = dir0 + 'raw/ParkerMacCready2017CTDDataFeb2018.xlsx'
        sheet_name = '2017Provisional_CTDResults'
    elif year == 2018:
        read_new = True
        #ctd_fn = dir0 + 'raw/Parker_2018.xlsx'
        ctd_fn = dir0 + 'raw/ParkerMacCready2018CTDDOMar2020.xlsx'
        sheet_name = '2018_CTDDOResults'
    elif year == 2019:
        read_new = True
        ctd_fn = dir0 + 'raw/ParkerMacCready2019CTDDataFeb2020.xlsx'
        sheet_name = '2019Provisional_CTDResults'
    else:
        ctd_fn = dir0 + 'raw/ParkerMacCready1999-2016CTDDataMay2018.xlsx'
        sheet_name = '1999-2016Finalized_CTDResults'

    # data long names; we retain only these fields
    
    # DEFAULTS
    date_col_name = 'Date'
    station_col_name = 'Station'
    depth_col_name = 'Depth'
    data_new_names = ['Salinity', 'Temperature', 'Sigma', 'Chl', 'DO', 'Turb', 'Z']
    
    if year == 2017:
        data_original_names = ['Salinity', 'Temp', 'Density',
                           'Chla_adjusted', 'DO_raw',
                           'Turbidity', 'Z']
    elif year == 2018: # missing Chl
        #date_col_name = 'UTCDate'
        #station_col_name = 'SiteCode'
        #depth_col_name = 'ActualDepthDecimal'
        data_original_names = ['Salinity', 'Temp', 'Density', 'DO_adjusted', 'Z']
        data_new_names = ['Salinity', 'Temperature', 'Sigma', 'DO', 'Z']
        
    elif year == 2019: # missing Chl
        data_original_names = ['Salinity', 'Temp', 'Density', 'DO_raw',
                           'Turbidity', 'Z']
        data_new_names = ['Salinity', 'Temperature', 'Sigma', 'DO', 'Turb', 'Z']
    else:
        data_original_names = ['Salinity', 'Temp', 'Density',
                           'Chla_adjusted', 'DO_adjusted',
                           'Turbidity', 'Z']
    
    # read in the data (all stations, all casts)
    if read_new:
        all_casts = pd.read_excel(ctd_fn, sheet_name=sheet_name,
            parse_dates = [date_col_name])
        read_new = False
    
        
    varname_dict = dict(zip(data_original_names, data_new_names))
    
    # initialize a DataFrame to hold all cast data for a year
    Casts = pd.DataFrame()

    # loop over all stations
    for station in sta_df.index:
        print(' - processing: ' + station)
        casts = all_casts[all_casts[station_col_name] == station]
        casts = casts.set_index(date_col_name)
        casts['Z'] = -casts[depth_col_name] # and make a Z column
        casts = casts[data_original_names] # keep only selected columns
        casts = casts.rename(columns=varname_dict) # rename columns
    
        # identify a single cast by its DATE
        alldates = casts.index
        alldates = alldates[alldates.year == year]
        castdates = alldates.unique() # a short list of unique dates (1 per cast)

        # process all casts for this station
        for cdate in castdates:
            cast = casts[casts.index==cdate]
            # drop repeat values (aiming for the first of a depth pair)
            zdf = np.diff(cast['Z'])
            zdf = np.concatenate((np.array([1.,]),zdf))
            mask = zdf != 0
            cast = cast[mask]
            cast = cast.sort_values('Z', ascending=True)
            # now packed bottom to top (like HYCOM and ROMS)
            cast = cast.iloc[5:,:] # drop deepest values (sometimes noisy)
            # add some columns that are needed later
            cast['Station'] = station
            cast['Date'] = cdate
            # We add Date as a column becasue we drop the index (Date)
            # in the concat operation below.  The resulting DataFrame "Casts"
            # just has numbers for its index
            Casts = pd.concat((Casts, cast), ignore_index=True)
    
    if year == 2018:
        Casts['Turb'] = np.nan
        Casts['Chl'] = np.nan
        
    # save result
    Casts.to_pickle(dir0 + 'Casts_' + str(year) + '.p')
