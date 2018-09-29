#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 10:04:20 2018

Plots data from a Department of Ecology, for all CTD stations,
and associated model extractions created by
LiveOcean/extract/driver_get_cast.py.

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

# set to False to save pngs, True to see on screen
testing = True
if testing:
    save_fig = False
else:
    save_fig = True

Ldir['gtagex'] = 'cas4_v1_lo6biom'

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

# where to put plots
if save_fig==True:
    dir11 = Ldir['parent'] + 'ptools_output/ecology/'
    Lfun.make_dir(dir11)
    dir1 = dir11 + 'casts' + str(year) + '_'+ Ldir['gtagex'] + '/'
    Lfun.make_dir(dir1, clean=True)

# choose which data fields to plot by commenting out parts of this list
data_to_plot = ['Salinity', 'Temperature']
    
# data short names, units, and plot ranges
data_names =  ['Salinity','Temperature','Sigma', 'Chl', 'DO',   'Turb', 'Z']
data_units =  ['psu',     'deg C',      'kg/m3', 'ug/l', 'mg/l', '',    'm']
data_ranges = [(14,34),   (4,25),       (14,26), (0,40), (0,18), (0,4),(-125,0)]

# dictionaries to look up data attributes using names
data_unit_dict = dict(zip(data_names, data_units))
data_range_dict = dict(zip(data_names, data_ranges))

# lists and dictionaries for plotting different months

months = range(1,13) # a list of 1 to 12
    
month_color_dict = dict(zip(months,
    ['mediumblue', 'royalblue', 'cadetblue', 'aquamarine',
    'lightgreen', 'greenyellow', 'gold', 'orange',
    'lightsalmon', 'mediumorchid', 'slateblue', 'purple']))
month_name_dict = dict(zip(months,
    ['Jan','Feb','Mar','Apr','May','Jun',
     'Jul','Aug','Sep','Oct','Nov','Dec']))

# plotting setup
plt.close('all')
figsize = (12,7)

# trim the station list if desired
if testing:
    sta_list = [sta for sta in sta_df.index if 'DNA001' in sta]
else:
    sta_list = [sta for sta in sta_df.index]

for station in sta_list:
    
    # set up the plot axes (an array)
    NR = 2
    NC = len(data_to_plot)
    fig, axes = plt.subplots(nrows=NR, ncols=NC, sharey=True,
            figsize=figsize, squeeze=False)

    # loop over all casts at this station
    print(' - Working on: ' + station)           
    casts = Casts[Casts['Station'] == station]   
    casts = casts.set_index('Date')    
    # identify a single cast by its DATE
    alldates = casts.index
    castdates = alldates.unique() # a short list of unique dates (1 per cast)
    title = str(year) + ' || ' + station + ': ' + sta_df.loc[station,'Descrip']
    Max_z = -float(sta_df.loc[station, 'Max_Depth'])
    ylim = data_range_dict['Z']
    
    # plot the CTD cast data for this station
    for cdate in castdates:
        imo = cdate.month
        
        if imo in months:
            
            cast = casts[casts.index==cdate]
            
            # also get the model fields for this day
            date_string = cdate.strftime('%Y.%m.%d')
            dir0m = Ldir['parent'] + 'LiveOcean_output/casts/'
            fnm = dir0m + Ldir['gtagex'] + '/' + station + '_' + date_string + '.nc'
            try:
                ds = nc4.Dataset(fnm)
                m_dict = dict()
                m_dict['Salinity'] = ds['salt'][:].squeeze()
                m_dict['Temperature'] = ds['temp'][:].squeeze()
                z = ds['z_rho'][:].squeeze()
                z = z - z[-1] # reference to sea level
                m_dict['Z'] = z
                ds.close()
                plot_mod = True
            except OSError:
                plot_mod = False
                pass
                
            nc = 0
            for fld in data_to_plot:
                ax = axes[0,nc]
                cast.plot(x=fld, y = 'Z', style='-', grid=True,
                          color=month_color_dict[imo], ax=ax, legend=False,
                          xlim=data_range_dict[fld], ylim=ylim)
                ax_str = fld + ' (' + data_unit_dict[fld] + ')'
                ax.set_xlabel('')
                if nc==0:
                    ax.set_ylabel('Z (m)')
                ax.text(.95, .05, ax_str, fontsize=14, fontweight='bold',
                horizontalalignment='right', transform=ax.transAxes)
        
                if plot_mod:
                    axm = axes[1, nc]
                    axm.plot(m_dict[fld], m_dict['Z'], '-',
                            color=month_color_dict[imo])
                    axm.set_xlim(data_range_dict[fld])
                    axm.set_ylim(ylim)
                    axm.grid(True)
                    if nc==0:
                        axm.text(.05, .05, Ldir['gtagex'], fontsize=12,
                                 fontweight='bold', transform=axm.transAxes)
                nc += 1
        
    # add month labels with colors
    ax = axes[0, 0]
    for imo in months:
        ax.text(.05, 1 - imo/13.,
            month_name_dict[imo], color=month_color_dict[imo],
            fontsize=12,
            verticalalignment='center', transform=ax.transAxes)

    fig.suptitle(title, fontsize=16, fontweight='bold')

    if save_fig:
        plt.savefig(dir1 + station + '.png')
        plt.close()
    else:
        plt.show()
