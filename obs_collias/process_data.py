"""
Process data from Collias archive CTD+Bottle casts, and
save the results for future use.  The goal is to organize the data
as close as possible to the way it is organized in process_casts.py
in obs_ecy.

"""

import pandas as pd
import numpy as np
import pickle

testing = False
if testing:
    year_list = [1954]
else:
    year_list = range(1932, 1976) # data go from 1932 to 1975

# read in stored station info (links station names to locations)
dir0 = '../../ptools_data/collias/'
sta_df = pd.read_pickle(dir0 + 'sta_df.p')

# read in bottle data
bottle_fn = dir0 + 'raw/Historic_Results.xlsx'
df0 = pd.read_excel(bottle_fn)

# Note: in the Collias data all the chemical concentrations are given in mg/l.

# Make a dict to rename variables, using names that are consistent
# (except for units at this point) with the recent-era Ecology data
col_dict = {'station':'Station', 'date':'Date', 'time': 'Time (UTC)', 'depth5': 'Depth (m)',
       'temp':'Temp. (deg C)', 'salinity':'Salinity', 'density':'Density (kg m-3)',
       'DO': 'DO (mg L-1)',
       'SiOH4D': 'SiOH4 (mg L-1)', 'NO2D': 'NO2 (mg L-1)', 'NO3D': 'NO3 (mg L-1)',
       'OPD': 'OPD', 'alk':'Alkalinity'}
       
df = df0.rename(columns=col_dict)

df['Z (m)'] = -df['Depth (m)']
       
# keep only selected columns
df = df.reindex(columns=['Station', 'Date', 'Time (UTC)', 'Z (m)', 'Temp. (deg C)', 'Salinity',
    'DO (mg L-1)', 'NO3 (mg L-1)', 'NO2 (mg L-1)', 'SiOH4 (mg L-1)'])
# add some useful columns
df['Letters'] = [s[:3] for s in df['Station'].values]
df['Numbers'] = [s[-3:] for s in df['Station'].values]
df['Num0'] = [s[0] for s in df['Numbers'].values] # defines which basin we are in
# set index so we can loop over years
allbot = df.set_index('Date')


# loop over years
for year in year_list:
    
    # initialize a DataFrame to hold all data for a year
    Bottles = pd.DataFrame()
    
    ybot = allbot[allbot.index.year == year]

    # loop over stations
    for station in sta_df.index:
        sb = ybot[ybot['Station']==station]
        
        sb['Date'] = sb.index.copy()
        # We add Date column becasue we drop the index (Date) below.
        
        # Initialize DataFrame that holds all cleaned casts at this station.
        # Cleaned means no repeat depths
        Sb = pd.DataFrame()
        
        print_info = True
        ncast = 0
        # loop over dates (meaning individual casts)
        for dd in sb.index.unique(): # Note: the index is still Date at this point
            ncast += 1
            sbd = sb[sb.index==dd] # just data from this Date
            sbdz = sbd.set_index('Z (m)')
            # drop repeated values (duplicate depths)
            sbdzu = sbdz[~sbdz.index.duplicated()]
            sbdzu = sbdzu.sort_index() # sort deepest to shallowest
            if (len(sbdz) != len(sbdzu)) and print_info:
                print('%s %s dropped %d repeat bottles' % (station, str(dd), len(sbdz)-len(sbdzu)))
                
            Sb = pd.concat((Sb,sbdzu)) # this has "Z (m)" as its index
            
        if print_info and (len(sb>0)):
            print('*** %s There were %d casts at this station ***\n' % (station, ncast))

        Sb['Z (m)'] = Sb.index.copy() # save because we drop it in the concat below
        Bottles = pd.concat((Bottles, Sb), ignore_index=True, sort=False)
        # The final DataFrame "Bottles" just has numbers for its index.
        
    # save result
    Bottles.to_pickle(dir0 + 'Bottles_' + str(year) + '.p')
