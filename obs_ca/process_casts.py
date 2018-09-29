"""
Process data from Canadian CTD casts, and
save the results for future use.  I had to download the casts by hand
from the website, searching for stations 27 and 42 in 2017.

"""

import pandas as pd
import numpy as np
import pickle
from datetime import datetime
import matplotlib.pyplot as plt
plt.close('all')

import os

# where the data is, and where results will be stored
dir0 = '../../ptools_data/canada/'

# NOTE: list should only have files for a single year

fn_list = [item for item in os.listdir(dir0+'raw/') if '.ctd' in item]

#fn_list_alt = ['2017-01-0124.ctd', '2017-01-0118.ctd']

# initialize a DataFrame to hold all cast info
sta_df = pd.DataFrame(#index=fn_list_alt,
    columns=['Station', 'Desig', 'Descrip', 'Basin', 'Max_Depth',
        'Latitude', 'Longitude','Datetime','File'])
sta_df.index.name = 'Filename' # later we replace this

# initialize a DataFrame to hold all cast data for a year
Casts = pd.DataFrame()

sta_list = [] # used to ensure that sta_df only has one entry per station location

for fn in fn_list:
    dfn = dir0 + 'raw/' + fn
    
    print('--working on ' + fn)
    f = open(dfn,'r', errors='ignore')

    count = 0
    get_units = False

    vn_list = []
    unit_list = []
    for line in f:
        if '*END OF HEADER' in line:
            data_start = count
        
        if 'START TIME' in line:
            LS = line.split()
            ymd = LS[4]
            dt = datetime.strptime(ymd,'%Y/%m/%d')
        
        if 'STATION' in line:
            LS = line.split()
            station = 'SOG'+str(LS[2])
            
        if 'WATER DEPTH' in line:
            LS = line.split()
            water_depth = float(LS[3])
        
        sign_dict = {'N':1,'S':-1,'E':1,'W':-1}
        if 'LATITUDE' in line:
            LS = line.split()
            degs = LS[2]
            mins = LS[3]
            sign = LS[4]
            lat = sign_dict[sign] * (float(degs) + float(mins)/60)
            
        if 'LONGITUDE' in line:
            LS = line.split()
            degs = LS[2]
            mins = LS[3]
            sign = LS[4]
            lon = sign_dict[sign] * (float(degs) + float(mins)/60)
        
        if '$TABLE: CHANNELS' in line:
            get_units = True
        if ('$END' in line) and get_units:
            get_units = False
        if get_units:
            if "'deg C (ITS90)'" in line:
                line = line.replace("'deg C (ITS90)'",'degC')
            LS = line.split()
            if LS[0].isdigit():
                LS[1] = LS[1].split(':')[0]
                if LS[1] in vn_list:
                    LS[1] = LS[1]+'2'
                vn_list.append(LS[1])
                unit_list.append(LS[2])
        
        count += 1
    f.close()
    
    if station not in sta_list:
        # store metadata in sta_df
        sta_df.loc[fn,'Station'] = station
        sta_df.loc[fn,'Desig'] = 'X'
        sta_df.loc[fn,'Descrip'] = 'BLANK'
        sta_df.loc[fn,'Basin'] = 'Strait of Georgia'
        sta_df.loc[fn,'Max_Depth'] = water_depth
        sta_df.loc[fn,'Latitude'] = lat
        sta_df.loc[fn,'Longitude'] = lon
        sta_df.loc[fn,'Datetime'] = dt.strftime('%Y.%m.%d')
        sta_df.loc[fn,'File'] = fn
        sta_list.append(station)
        
    # get the cast data
    df = pd.read_csv(dfn,skiprows=data_start+1,
        header=None, names=vn_list, delim_whitespace=True)
    df.rename(columns={'Depth':'Z'}, inplace=True)
    df['Z'] = -df['Z']
    cast = df.copy()
    cast = cast.sort_values('Z', ascending=True)
    # now packed bottom to top (like HYCOM and ROMS)
    cast['Station'] = station
    cast['Date'] = dt
    # The resulting DataFrame "Casts" just has numbers for its index.
    Casts = pd.concat((Casts, cast), ignore_index=True)
    
    if False:
        # plotting
        z = cast['Z'].values
        temp = cast['Temperature'].values
        salt = cast['Salinity'].values
        fig, axes = plt.subplots(1,2,figsize=(12,6), sharey=True, squeeze=False)
        ax=axes[0,0]
        ax.plot(salt,z)
        ax.set_xlabel('Salinity')
        ax.set_ylabel('Z (m)')
        ax=axes[0,1]
        ax.plot(temp,z)
        ax.set_xlabel('Temperature (deg C)')
        fig.suptitle(fn)
        plt.show()
        
# reorganizing
sta_df = sta_df.set_index('Station')

# save result
year = dt.year
Casts.to_pickle(dir0 + 'Casts_' + str(year) + '.p')
sta_df.to_pickle(dir0 + 'sta_df.p')