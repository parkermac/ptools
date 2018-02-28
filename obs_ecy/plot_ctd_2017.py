#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 10:04:20 2018

@author: pm7
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 15:44:00 2018

Plots data from a Department of Ecology csv file, for flight CTD stations.

Designed to work with the new  2017 data I asked for.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
#import netCDF4 as nc4

## USER INPUT ##

# where are the data csv files
dir0 = '/Users/pm7/Documents/ptools_data/ecology/'

sta_df = pd.read_pickle(dir0 + 'sta_df.p')

# choose which station to plot
sta = 'ParkerMacCready2017CTDDataFeb2018'

fn = dir0 + sta + '.xlsx'

# choose which data fields to plot by commenting out parts of this list
data_to_plot = [
    'Salinity',
    'Temperature',
#    'Sigma',
#    'Chl',
#    'DO',
#    'Turb',
    ]

year = 2017 # integer year, or None to plot all years

## END USER INPUT ##

# lists of data properties

# data long names
# and we retain only these fields
data_long_names = ['Salinity', 'Temp', 'Density',
                   'Chla_adjusted', 'DO_raw',
                   'Turbidity', 'Z']

# data short names, units, and plot ranges
data_names =  ['Salinity','Temperature','Sigma', 'Chl', 'DO',   'Turb', 'Z']
data_units =  ['psu',     'deg C',      'kg/m3', 'ug/l', 'mg/l', '',    'm']
data_ranges = [(14,34),   (4,25),       (14,26), (0,40), (0,18), (0,3),(-100,0)]

# dictionaries to look up data attributes using names
data_name_dict = dict(zip(data_long_names, data_names))
data_unit_dict = dict(zip(data_names, data_units))
data_range_dict = dict(zip(data_names, data_ranges))

# more useful lists and dictionaries for plotting different months

months = range(1,13) # a list of 1 to 12

month_color_dict = dict(zip(months,
    ['mediumblue', 'royalblue', 'cadetblue', 'aquamarine',
    'lightgreen', 'greenyellow', 'gold', 'orange',
    'lightsalmon', 'mediumorchid', 'slateblue', 'purple']))

month_name_dict = dict(zip(months,
    ['Jan','Feb','Mar','Apr','May','Jun',
     'Jul','Aug','Sep','Oct','Nov','Dec']))


all_casts = pd.read_excel(fn, sheet_name='2017Provisional_CTDResults', parse_dates = ['Date'])

# setup
plt.close('all')
figsize = (18,7)

for station in sta_df.index:
    
    casts = all_casts[all_casts['Station'] == station]
    
    casts = casts.set_index('Date')

    casts['Z'] = -casts['Depth'] # and make a Z column
    casts = casts[data_long_names] # keep only selected columns
    casts = casts.rename(columns=data_name_dict) # and rename them
    
    # plot CASTS
    
    # identify a single cast by its DATE
    alldates = casts.index
    castdates = alldates.unique() # a short list of unique dates (1 per cast)
    
    title = sta
    if year != None:
        title = sta + ': ' + str(year)
        dt0 = datetime(year,1,1)
        dt1 = datetime(year,12,31)
        castdates = castdates[castdates >= dt0]
        castdates = castdates[castdates <= dt1] 
        casts = casts.loc[dt0:dt1,:]
        
    # set up the plot axes (an array)
    NR = 2
    NC = len(castdates)
    
    fig, axes = plt.subplots(nrows=NR, ncols=NC, sharey=True,
                             figsize=figsize, squeeze=False)    
    for cd in castdates:
        print('\'' + datetime.strftime(cd, '%Y.%m.%d') + '\',')
            
    # plot the CTD cast data for this station
    for cdate in castdates:
        imo = cdate.month
        cast = casts[casts.index==cdate]
        # drop repeat values (aiming for the first of a depth pair)
        zdf = np.diff(cast['Z'])
        zdf = np.concatenate((np.array([1.,]),zdf))
        mask = zdf != 0
        cast = cast[mask]
        cast = cast[:-5] # drop bottom values (sometimes bad)
            
        for fld in data_to_plot:
            if fld == 'Salinity':
                nr = 0
            elif fld == 'Temperature':
                nr = 1
                
            nc = imo-1
            
            ax = axes[nr,nc]
            cast.plot(x=fld, y = 'Z', style='--', grid=True,
                      color=month_color_dict[imo], ax=ax, legend=False,
                      xlim=data_range_dict[fld], ylim=data_range_dict['Z'])
            
            #ax.plot(m_dict[fld], m_dict['Z'], '-', color=month_color_dict[imo])
            
            ax.set_xlabel(fld + ' (' + data_unit_dict[fld] + ')')
            
    
    # set more things about the cast plots
    
    # add month labels with colors
    if nr == 0:
        ax = axes[0, 0]
        for imo in months:
            ax.text(.05, 1 - imo/13.,
                month_name_dict[imo], color=month_color_dict[imo],
                fontsize=14, fontweight='bold',
                verticalalignment='center', transform=ax.transAxes)
    
    fig.suptitle(title)

    
    plt.show()
