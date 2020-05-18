#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots data from the Department of Ecology, combining data from CTD casts
and bottles at the same station.  Also including model output.  4 depths.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import netCDF4 as nc4

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart()
import zfun

sys.path.append(os.path.abspath(Ldir['LO'] + 'plotting'))
import pfun

# ***** User Edits

# specify which model run to use
Ldir['gtagex'] = 'cas6_v3_lo8b'

testing = False

do_map = True

for_web = True

year_list = [2017, 2018, 2019]

# +++ load ecology CTD info +++
dir0 = Ldir['parent'] + 'ptools_data/ecology/'
sta_df = pd.read_pickle(dir0 + 'sta_df.p')
try:
    # add Canadian info
    dir1 = Ldir['parent'] + 'ptools_data/canada/'
    sta_df_ca = pd.read_pickle(dir1 + 'sta_df.p')
    sta_df = pd.concat((sta_df, sta_df_ca), sort=False)
except FileNotFoundError:
    pass

fs=16 # font size

if testing:
    year_list = [2017]
    sta_to_plot = [s for s in sta_df.index if 'SAR003' in s]
    save_fig = False
else:
    sta_to_plot = [s for s in sta_df.index]
    save_fig = True
    
if for_web:
    fs = 10
    tag = 'web'
    do_map = False
    figsize = (6,4)
    NR = 2; NC = 2
else:
    tag = 'val'
    if do_map:
        figsize = (20,8)
        NR = 2; NC = 3
    else:
        figsize = (15,8)
        NR = 2; NC = 2

print('do_map = ' + str(do_map))
print('for_web = ' + str(for_web))

# ***** End User Edits

# where to put plots
dir11 = Ldir['parent'] + 'ptools_output/ecology/'

# # plotting
plt.rc('font', size=fs)
plt.close('all')

# y limits for plotting
lim_dict = {'Salinity': (10,35), 'Temp. (deg C)': (5,20),
'DO (mg L-1)': (0, 15), 'DIN (uM)': (0,40)}

# lim_dict = {'Salinity': (0, 34), 'Temp. (deg C)': (0, 24),
#     'DO (mg L-1)': (0, 20), 'DIN (uM)': (0,55), 'Chl (mg m-3)': (0,50)}
lab_dict = {'Salinity':'(a) Salinity', 'Temp. (deg C)':'(b) Temperature [$^{\circ}C$]',
    'DO (mg L-1)': '(c) DO [$mg \ L^{-1}$]', 'DIN (uM)':'(d) DIN  [$\mu M$]'}

# line colors (by depth)
clist = ['r', 'g', 'b', 'k']

for year in year_list:
    
    print('\n'+str(year)+30*'-=')
    
    if save_fig==True:
        if do_map == False:
            dir1 = dir11 + tag + '_series_' + Ldir['gtagex'] + '_'+ str(year) + '_nomap/'
        else:
            dir1 = dir11 + tag + '_series_' + Ldir['gtagex'] + '_'+ str(year) + '/'
        Lfun.make_dir(dir1, clean=True)
    
    in_fn = 'ObsMod_' + Ldir['gtagex'] + '_'+ str(year) + '.p'
    print('Loading ' + in_fn)
    Bc = pd.read_pickle(dir11 + in_fn)
    Z_list = list(Bc['Znom'].unique())
    
    # x limits for plotting
    dt0 = pd.datetime(year,1,1)
    dt1 = pd.datetime(year,12,31)
    
    # loop over all stations
    for station in sta_to_plot:
        fig = plt.figure(figsize=figsize)
        A = Bc[Bc['Station'] == station]   
        A = A.set_index('Date')
        pp = 1
        
        for vn in ['Salinity', 'Temp. (deg C)', 'DO (mg L-1)', 'DIN (uM)']:
            ax = fig.add_subplot(NR, NC ,pp)
            ii = 0
            for Zn in Z_list:
                a = A[A['Znom']==Zn]
                try:
                    a.plot(y=[vn, 'Mod '+vn], ax=ax,
                        style=['-o'+clist[ii], '--+'+clist[ii]], legend=False)
                except:
                    pass
                ii += 1
            ax.set_xlim(dt0,dt1)
            ax.set_ylim(lim_dict[vn])
            ax.xaxis.set_major_locator(mdates.MonthLocator())
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
            ax.xaxis.set_tick_params(labelrotation=45)
            ax.text(.05, .9, lab_dict[vn], fontweight='bold', transform=ax.transAxes)
            if pp==1:
                ax.text(.05, .05, 'Solid = Observed, Dashed = Model',
                    style='italic', color='k', transform=ax.transAxes)
            if pp == 3:
                z_counter = 0
                for zz in Z_list:
                    ax.text(.05, .3 - .08*z_counter, 'Z =   %d m' % (zz),
                        color=clist[z_counter], fontweight='bold', transform=ax.transAxes)
                    z_counter += 1
            if pp in [1,2]:
                ax.set_xticklabels('')
                ax.set_xlabel('')
            else:
                ax.set_xlabel('Date ' + str(year))
            ax.grid(True)
            pp += 1
            if (pp==3) and do_map:
                pp = 4
    
        if do_map:
            # add station location map
            ax = fig.add_subplot(1,3,3)
            lon = sta_df.loc[station, 'Longitude']
            lat = sta_df.loc[station, 'Latitude']
            ax.plot(lon, lat, '*r')
            ax.text(lon+.01, lat, station, color='b', fontsize=12, fontweight='bold')
            ax.set_title(sta_df.loc[station, 'Descrip'])
            pfun.add_coast(ax)
            pfun.dar(ax)
            if sta_df.loc[station, 'Basin'] in ['Grays Harbor', 'Willapa Bay']:
                # Coastal Estuaries
                ax.set_xlim(-124.3, -123.6)
                ax.set_ylim(46.3, 47.1)
            else:
                # Puget Sound
                ax.set_xlim(-124, -122)
                ax.set_ylim(47, 49.5)
        
        fig.tight_layout()
    
        if save_fig:
            plt.savefig(dir1 + station + '.png')
            plt.close()
        else:
            plt.show()

plt.rcdefaults()
