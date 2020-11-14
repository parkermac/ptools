"""
Plot one or more CTD casts with an emphasis on stratification.

For Julie K's Puget Sound class, 2020_10

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart()

sys.path.append(os.path.abspath('../../LiveOcean/plotting'))
import pfun

# +++ load ecology CTD cast data +++
dir0 = Ldir['parent'] + 'ptools_data/ecology/'
# load processed station info and data
sta_df = pd.read_pickle(dir0 + 'sta_df.p')
year = 2017
Casts = pd.read_pickle(dir0 + 'Casts_' + str(year) + '.p')


# choose which data fields to plot by commenting out parts of this list
data_to_plot = ['Salinity', 'Temperature', 'Sigma']
    
# data short names, units, and plot ranges
data_names =  ['Salinity','Temperature','Sigma', 'Chl', 'DO',   'Turb', 'Z']
data_units =  ['psu',     'deg C',      'kg/m3', 'ug/l', 'mg/l', '',    'm']
data_ranges = [(14,34),   (4,25),       (14,26), (0,40), (0,18), (0,4),(-125,0)]

# dictionaries to look up data attributes using names
data_unit_dict = dict(zip(data_names, data_units))
data_range_dict = dict(zip(data_names, data_ranges))

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
fs = 14
plt.rc('font', size=fs)
lw = 3
aa = [-124, -122, 47, 49]

outdir = Ldir['parent']+'ptools_output/ecology/stratification/'
Lfun.make_dir(outdir, clean=True)

testing = False
if testing:
    sta_list = ['PSB003']
else:
    sta_list = list(sta_df.index)
    
for station in sta_list:

    fig = plt.figure(figsize=(20,8))
    ax1 = fig.add_subplot(141)
    ax2 = fig.add_subplot(142)
    ax3 = fig.add_subplot(143)
    ax4 = fig.add_subplot(144) # map

    casts = Casts[Casts['Station'] == station]   
    casts = casts.set_index('Date')    
    # identify a single cast by its DATE
    alldates = casts.index
    castdates = alldates.unique() # a short list of unique dates (1 per cast)
    Max_z = -float(sta_df.loc[station, 'Max_Depth'])
    # plot the CTD cast data for this station
    for cdate in castdates:
        imo = cdate.month
        cast = casts[casts.index==cdate]
        cast.plot(x='Salinity', y='Z', ax=ax1, legend=False, grid=True, lw=lw, color=month_color_dict[imo])
        cast.plot(x='Temperature', y='Z', ax=ax2, legend=False, grid=True, lw=lw, color=month_color_dict[imo])
        cast.plot(x='Sigma', y='Z', ax=ax3, legend=False, grid=True, lw=lw, color=month_color_dict[imo])
    
    ax2.set_yticklabels([])
    ax3.set_yticklabels([])
    ax1.set_ylabel(r'Z $[m]$')

    ax1.set_xlabel(r'Salinity $[g\ kg^{-1}]$')
    ax2.set_xlabel(r'Temperature $[^{\circ} C]$')
    ax3.set_xlabel(r'Density - 1000 $[kg\ m^{-3}]$')


    ax1.set_ylim(top=0)
    ax2.set_ylim(top=0)
    ax3.set_ylim(top=0)

    for imo in months:
        ax1.text(.05, 1 - imo/13.,
            month_name_dict[imo], color=month_color_dict[imo],
            verticalalignment='center', transform=ax1.transAxes, weight='bold')

    fig.suptitle(station + ': ' + sta_df.loc[station, 'Descrip'], weight='bold')
    
    for sn in sta_df.index:
        lon = sta_df.loc[sn, 'Longitude']
        lat = sta_df.loc[sn, 'Latitude']
        if lon > aa[0] and lat > aa[2]:
            if sn == station:
                ax4.plot(lon, lat, '*r')
                ax4.text(lon+.01, lat, sn, color='r', weight='bold')
            else:
                ax4.plot(lon, lat, '*b', alpha=.4)
                ax4.text(lon+.01, lat, sn, size=6, color='b', weight='bold', alpha=.4)
    # Puget Sound
    ax4.axis(aa)
    ax4.set_xticks([-124, -123, -122])
    ax4.set_yticks([47, 48, 49])
    pfun.add_coast(ax4)
    pfun.dar(ax4)
    
    if testing:
        plt.show()
    else:
        plt.savefig(outdir + station + '.png')
        plt.close()
    
plt.rcdefaults()