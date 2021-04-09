"""
This reads in and processes the WOAC csv files from Simone Alin,
specifically to create a sta_df.p file that has all the station
info needed for a cast extraction.

"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun

in_dir = Path(__file__).absolute().parent.parent.parent / 'ptools_data' / 'woac' / 'raw_2019_10'

out_dir = Path(__file__).absolute().parent.parent.parent / 'ptools_output' / 'woac2'
Lfun.make_dir(out_dir)

cruise_list = ['AQ201710', 'CAB1079', 'NORSEMANIIOCT2018', 'RBTSN201805', 'RC001', 'RC006', 'RC007']

# cruise_list = ['RC001', 'RC007'] # for testing

for cruise in cruise_list:
    print(cruise.center(60,'-'))
    
    # read the file to a DataFrame and combine two columns to start the datetime
    # as a column called DATE_UTC_TIME_UTC
    fn = in_dir / (cruise + '_data.csv')
    df = pd.read_csv(fn, parse_dates=[['DATE_UTC','TIME_UTC']])
    # the index at this point is just row number
    
    # set known bad data to np.nan
    df[df=='nan nan'] = np.nan
    df[df==-999] = np.nan
    
    # drop rows or columns that have no good data at all, or that are empty
    df = df.dropna(how='all')
    
    # keep only selected columns
    df = df[['LONGITUDE_DEC', 'LATITUDE_DEC', 'STATION_NO','DATE_UTC_TIME_UTC']]
    # and rename them for convenience
    df = df.rename(columns={'LONGITUDE_DEC':'Longitude', 'LATITUDE_DEC':'Latitude',
        'STATION_NO': 'Station', 'DATE_UTC_TIME_UTC':'Datetime'})
    
    # convert the datetime to a float, so that is it not lost in the groupby mean()
    df['Datetime'] = pd.to_datetime(df['Datetime'])
    df['Datetime'] = df['Datetime'].values.astype(float)
    
    # check that a station number was only visited once in a cruise
    for sta in df['Station'].unique():
        a = df[df['Station']==sta]
        b = np.diff(a.index)
        if not (b==1).all():
            print('** Warning: repeat Station %d in Cruise %s' % (sta, cruise))
    
    # take mean of everything by station number
    dfg = df.groupby('Station').mean()
    # now the index is station number
    
    # convert back to datetime, and round to the nearest hour
    dfg['Datetime'] = pd.to_datetime(dfg['Datetime']).dt.round('H')
    
    # make sure that sign of Longitude is negative
    dfg['Longitude'] = -np.abs(dfg['Longitude'])
    
    # add the cuise name as a column
    dfg['Cruise'] = cruise
    
    if cruise == cruise_list[0]:
        A = dfg.copy()
    else:
        A = pd.concat((A, dfg))
        
# save output
out_fn = out_dir / 'sta_df.p'
print('Saving to:\n' + str(out_fn))
A.to_pickle(out_fn)

# plot station locations
plt.close('all')
A.plot(x='Longitude', y = 'Latitude', marker='*', linestyle='')
plt.show()
    
