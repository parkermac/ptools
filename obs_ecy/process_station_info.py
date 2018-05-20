"""
Process station info from WA Dept. of Ecology CTD casts, and
save the results for future use.

"""

import pandas as pd

# where the data is, and where results will be stored
dir0 = '../../ptools_data/ecology/'

# file with station info (good for all years?)
sta_fn = dir0 + 'raw/ParkerMacCreadyCoreStationInfoFeb2018.xlsx'

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
sta_df.to_pickle(dir0 + 'sta_df.p')

