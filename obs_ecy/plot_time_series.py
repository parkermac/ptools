#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots data from a Department of Ecology station, and
a model mooring extraction, as time series.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()

testing = True
Ldir['gtagex'] = 'cas4_v0_lo6m'

# +++ load ecology CTD cast data +++
dir0 = Ldir['parent'] + 'ptools_data/ecology/'
# load processed station info and data
sta_df = pd.read_pickle(dir0 + 'sta_df.p')
# add Canadian data
dir1 = Ldir['parent'] + 'ptools_data/canada/'
# load processed station info and data
sta_df_ca = pd.read_pickle(dir1 + 'sta_df.p')
sta_df = pd.concat((sta_df, sta_df_ca))
year = 2017
Casts = pd.read_pickle(dir0 + 'Casts_' + str(year) + '.p')
Casts_ca = pd.read_pickle(dir1 + 'Casts_' + str(year) + '.p')
Casts = pd.concat((Casts, Casts_ca))

if testing==True:
    sta_to_plot = [s for s in sta_df.index if 'PSB003' in s]
    save_fig = False
else:
    sta_to_plot = [s for s in sta_df.index]
    save_fig = True

# where to put plots
if save_fig==True:
    dir11 = Ldir['parent'] + 'ptools_output/ecology/'
    Lfun.make_dir(dir11)
    dir1 = dir11 + 'series' + str(year) + '_' + Ldir['gtagex'] + '/'
    Lfun.make_dir(dir1, clean=True)
        
# initialize dicts to hold all time series
temp_dict = dict()
salt_dict = dict()

for station in sta_to_plot:
    
    # loop over all casts at this station
    print(' - Working on: ' + station)           
    casts = Casts[Casts['Station'] == station]   
    casts = casts.set_index('Date')    
    #casts = casts.rename(columns=data_name_dict) # rename columns
    # identify a single cast by its DATE
    alldates = casts.index
    castdates = alldates.unique() # a short list of unique dates (1 per cast)
    
    # initialize DataFrames for depth-range-averaged time series
    # at this station
    temp_df = pd.DataFrame(index=castdates, columns=['sm','dm','so','do'])
    salt_df = pd.DataFrame(index=castdates, columns=['sm','dm','so','do'])
    # NAMING CONVENTION
    # so = shallow obs
    # do = deep obs
    # sm = shallow model
    # dm = deep model
    
    # gather data for time series
    for cdate in castdates:
        imo = cdate.month
        cast = casts[casts.index==cdate]
        
        temp_df.loc[cdate, 'so'] = cast['Temperature'][cast['Z']>-5].mean()
        temp_df.loc[cdate, 'do'] = cast['Temperature'][cast['Z']<=-5].mean()
        salt_df.loc[cdate, 'so'] = cast['Salinity'][cast['Z']>-5].mean()
        salt_df.loc[cdate, 'do'] = cast['Salinity'][cast['Z']<=-5].mean()
    
        # also get the model fields for this day
        date_string = cdate.strftime('%Y.%m.%d')
        dir0m = Ldir['LOo'] + 'casts/'
        fnm = dir0m + Ldir['gtagex'] + '/' + station + '_' + date_string + '.nc'
        try:
            ds = nc4.Dataset(fnm)
            salt = ds['salt'][:].squeeze()
            temp = ds['temp'][:].squeeze()
            z = ds['z_rho'][:].squeeze()
            z = z - z[-1] # reference to sea level
            ds.close()
            temp_df.loc[cdate, 'sm'] = temp[z>-5].mean()
            temp_df.loc[cdate, 'dm'] = temp[z<=-5].mean()
            salt_df.loc[cdate, 'sm'] = salt[z>-5].mean()
            salt_df.loc[cdate, 'dm'] = salt[z<=-5].mean()
        except OSError:
            pass
            
    temp_dict[station] = temp_df
    salt_dict[station] = salt_df
    
#PLOTTING
plt.close('all')

dt0 = pd.datetime(2017,1,1)
dt1 = pd.datetime(2017,12,31)
col_vec = ['so', 'sm', 'do', 'dm']
name_vec = ['Shallow Obs', 'Shallow Model', 'Deep Obs', 'Deep Model']
style_vec = ['-or', '--+r', '-ob', '--+b']
name_dict = dict(zip(col_vec, name_vec))
style_dict = dict(zip(col_vec, style_vec))

lw=2
for station in sta_to_plot:
    Temp_df = temp_dict[station]
    Salt_df = salt_dict[station]
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12,7), squeeze=False)
    fig.suptitle(station + ' ' + Ldir['gtagex'], fontsize=16, fontweight='bold')
    
    ax = axes[0,0]
    for cc in col_vec:
        try:
            Temp_df.plot(y=cc, style=style_dict[cc], label=name_dict[cc], linewidth=lw, ax=ax)
        except TypeError:
            pass
    
    ax.set_xlim(dt0,dt1)
    ax.set_ylim(4,20)
    # ax.set_xticklabels('')
    # ax.set_xlabel('')
    ax.grid(True)
    ax.text(.5,.05,'Temperature (deg C)', horizontalalignment='center',
        fontweight='bold', transform=ax.transAxes)
    
    ax = axes[0,1]
    for cc in col_vec:
        try:
            Salt_df.plot(y=cc, style=style_dict[cc], linewidth=lw, ax=ax, legend=False)
        except TypeError:
            pass
    #ax.set_xlabel('Date')
    ax.set_xlim(dt0,dt1)
    ax.set_ylim(8,34)
    ax.text(.05,.1,'Salinity', fontweight='bold', transform=ax.transAxes)
    ax.grid(True)
    
    if save_fig:
        plt.savefig(dir1 + station + '.png')
        plt.close()
    else:
        plt.show()
    
