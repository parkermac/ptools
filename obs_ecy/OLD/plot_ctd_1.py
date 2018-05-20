#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 15:44:00 2018

Plots data from a Department of Ecology csv file, for flight CTD stations.

Designed to work with the new HCB010 2017 data I asked for.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import netCDF4 as nc4

## USER INPUT ##

# where are the data csv files
dir0 = '/Users/pm7/Documents/ptools_data/ecology/'

# choose which station to plot
sta = 'HCB010_2017_custom'

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


casts = pd.read_excel(fn, sheet_name='Data', parse_dates = ['Date'])
casts = casts.set_index('Date')
casts['Z'] = -casts['Depth'] # and make a Z column
casts = casts[data_long_names] # keep only selected columns
casts = casts.rename(columns=data_name_dict) # and rename them

# PLOTTING

# setup
plt.close('all')
figsize = (18,7)

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
    
# intitialize data frames to save time series data
ser_df_top = pd.DataFrame(index=castdates, columns=data_to_plot)
ser_df_bot = pd.DataFrame(index=castdates, columns=data_to_plot)

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
    
    # also get the model fields for this day
    date_string = datetime.strftime(cdate, '%Y.%m.%d')
    dir0m = '/Users/pm7/Documents/LiveOcean_output/cast/'
    fnm = dir0m + 'cas3_v0_lo6m_HCB010_' + date_string + '.nc'
    ds = nc4.Dataset(fnm)
    m_dict = dict()
    m_dict['Salinity'] = ds['salt'][:].squeeze()
    m_dict['Temperature'] = ds['temp'][:].squeeze()
    m_dict['Z'] = ds['z_rho'][:].squeeze()
    ds.close()

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
        
        ax.plot(m_dict[fld], m_dict['Z'], '-', color=month_color_dict[imo])
        
        ax.set_xlabel(fld + ' (' + data_unit_dict[fld] + ')')
        
        # also gather a time series entry
        zdiv = -10
        #ser_df_top.loc[cdate, fld] = cast[(cast['Z']>=zdiv) & (cast['Z']<-3)].mean(axis=0)[fld]
        ser_df_top.loc[cdate, fld] = cast[cast['Z']>=zdiv].mean(axis=0)[fld]
        ser_df_bot.loc[cdate, fld] = cast[cast['Z']<zdiv].mean(axis=0)[fld]

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

if False:
    # NEXT: A figure for TIME SERIES
    
    # make the axes (overwrites previous axes object)
    fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=figsize,
                             squeeze=False, sharex=True)
    
    # plot data
    nr = 0
    dft = ser_df_top
    dfb = ser_df_bot
    for fld in data_to_plot:
        nc = get_nc[fld] 
        ax = axes[nr,nc]
        dft.plot(ax=ax, y=fld, grid=True, style='*-r', legend=False)
        dfb.plot(ax=ax, y=fld, grid=True, style='*-b', legend=False,
                 ylim=data_range_dict[fld])
        ax.text(.05, .95, fld + ' (' + data_unit_dict[fld] + ')',
                fontsize=12, fontweight='bold', transform=ax.transAxes)
        ax.set_xlabel('')
    
    # add explanatory text
    ax = axes[0, 0]
    ax.text(.05, .1, 'Mean above ' + str(zdiv) + ' m', color='r',
        fontsize=14, transform=ax.transAxes)
    ax.text(.05, .05, 'Mean below ' + str(zdiv) + ' m', color='b',
        fontsize=14, transform=ax.transAxes)
    
    fig.suptitle(title)

plt.show()
