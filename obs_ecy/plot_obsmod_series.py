#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots data from the Department of Ecology, combining data from CTD casts
and bottles at the same station.  Also including model output.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import netCDF4 as nc4

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()
import zfun

pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

# ***** User Edits

year = 2017

# specify which model run to use
#Ldir['gtagex'] = 'cas4_v2_lo6biom'
#Ldir['gtagex'] = 'cas5_v3_lo8'
Ldir['gtagex'] = 'cas6_v3_lo8b'

testing = False
no_map = True

for_web = False
if for_web:
    no_map = True

# ***** End User Edits

# +++ load ecology CTD cast data +++
dir0 = Ldir['parent'] + 'ptools_data/ecology/'
# load processed station info and data
sta_df = pd.read_pickle(dir0 + 'sta_df.p')
Casts = pd.read_pickle(dir0 + 'Casts_' + str(year) + '.p')
try:
    # add Canadian data
    dir1 = Ldir['parent'] + 'ptools_data/canada/'
    # load processed station info and data
    sta_df_ca = pd.read_pickle(dir1 + 'sta_df.p')
    sta_df = pd.concat((sta_df, sta_df_ca), sort=False)
    Casts_ca = pd.read_pickle(dir1 + 'Casts_' + str(year) + '.p')
    Casts = pd.concat((Casts, Casts_ca), sort=False)
except FileNotFoundError:
    pass

# +++ load ecology bottle data +++
try:
    Bottles = pd.read_pickle(dir0 + 'Bottles_' + str(year) + '.p')
    do_bottles = True
except FileNotFoundError:
    do_bottles = False

# still need to get Canadian bottle data

if testing==True:
    sta_to_plot = [s for s in sta_df.index if 'SAR003' in s]
    save_fig = False
else:
    sta_to_plot = [s for s in sta_df.index]
    save_fig = True
    if for_web==True:
        tag = 'web'
    else:
        tag = 'val'

# where to put plots
if save_fig==True:
    dir11 = Ldir['parent'] + 'ptools_output/ecology/'
    Lfun.make_dir(dir11)
    if no_map == True:
        dir1 = dir11 + tag + '_series_' + Ldir['gtagex'] + '_'+ str(year) + '_nomap/'
    else:
        dir1 = dir11 + tag + '_series_' + Ldir['gtagex'] + '_'+ str(year) + '/'
    Lfun.make_dir(dir1, clean=True)

Bc = pd.DataFrame() # DataFrame to hold all bottle and cast data

# loop over all stations
for station in sta_to_plot:
    
    casts = Casts[Casts['Station'] == station]
    casts = casts.set_index('Date')    
    # identify a single cast by its DATE
    calldates = casts.index
    castdates = calldates.unique() # a short list of unique dates (1 per cast)
    # also for bottles
    if do_bottles:
        bottles = Bottles[Bottles['Station'] == station]   
        bottles = bottles.set_index('Date')    
        # identify a single cast by its DATE
        balldates = bottles.index
        bottledates = balldates.unique() # a short list of unique dates (1 per cast)
        # all valid dates (really bottledates should be a subset of castdates, right?)
        dates = bottledates.union(castdates)
    else:
        dates = castdates
    
    # some useful dictionaries for renaming and rescaling
    if year == 2019:
        cast_vn_dict = {'Salinity': 'Salinity', 'Temperature': 'Temp. (deg C)',
            'DO': 'DO (mg L-1)'}
    else:
        cast_vn_dict = {'Salinity': 'Salinity', 'Temperature': 'Temp. (deg C)',
            'Chl': 'Chl (mg m-3)', 'DO': 'DO (mg L-1)'}
    bot_vn_dict = {'DIN': 'DIN (uM)'}
    mod_fac_dict = {'salt': 1, 'temp': 1,
        'NO3': 1, 'phytoplankton': 2.5, 'oxygen': 0.032}
    # model fields have these units AFTER multiplication by mod_fac_dict
    mod_vn_dict = {'salt': 'Mod Salinity', 'temp': 'Mod Temp. (deg C)',
        'NO3': 'Mod DIN (uM)', 'phytoplankton': 'Mod Chl (mg m-3)', 'oxygen': 'Mod DO (mg L-1)'}
    
    # rename the columns
    casts = casts.rename(columns=cast_vn_dict)
    if do_bottles:
        bottles = bottles.rename(columns=bot_vn_dict)
    
    # limit the columns
    casts = casts[['Station', 'Z'] + list(cast_vn_dict.values())]
    if do_bottles:
        bottles = bottles[['Station', 'Z', 'Znom'] + list(bot_vn_dict.values())]
    
    # Set the depths to look at.
    # [0,-10,-30] are the bottle nominal depths, but you can get deeper cast data
    # by adding to this list, e.g. -100
    Z_list = [0,-10,-30, -80]

    # loop over all casts at this station
    dates = dates[dates.year==year]
    
    for dd in dates:

        #Initialize a DataFrame for this station/date
        bc_columns = ['Station', 'Date', 'Znom', 'Z',
            'Salinity', 'Temp. (deg C)','Chl (mg m-3)', 'DO (mg L-1)','DIN (uM)',
            'Mod Salinity', 'Mod Temp. (deg C)', 'Mod Chl (mg m-3)',
            'Mod DO (mg L-1)', 'Mod DIN (uM)']
        bc = pd.DataFrame(index=Z_list, columns=bc_columns)
        bc['Date'] = dd
        bc['Station'] = station
        bc['Znom'] = Z_list # we save this as a column and as the index, because
                # later we will need it after we drop the index during
                # concatenation into Bc.
    
        # NOTE the brackets around [dd] keep the result as a DataFrame even if
        # we are only pulling out a single row.
        ca = casts.loc[[dd],:]
        ca = ca.set_index('Z')
        for Zn in Z_list:
            try:
                if Zn==0:
                    # OK to extrapolate at surface because shallowest cast data is
                    # always deeper than 0.
                    exn = False
                else:
                    # but don't extrapolate beyond the deepest cast data when looking
                    # for a given Zn
                    exn = True
                i0, i1, fr = zfun.get_interpolant(np.array([Zn]), ca.index.values, extrap_nan=exn)
                for vn in cast_vn_dict.values():
                    cv = ca[vn]
                    bc.loc[Zn,vn] = (1-fr)*cv.iloc[int(i0)] + fr*cv.iloc[int(i1)]
            except:
                pass
    
        if do_bottles:
            try:
                bo = bottles.loc[[dd],:]
                bo = bo.set_index('Znom')
                # for bottles we just put the data at Znom, ignoring the fact that
                # sometimes it is at a different actual depth.  When I looked at the
                # data it did not seem like a big problem.
                for Zn in Z_list:
                    for vn in bot_vn_dict.values():
                        bv = bo.loc[Zn, vn]
                        bc.loc[Zn,vn] = bv
            except KeyError:
                pass

        # finally get the model fields for this day
        date_string = dd.strftime('%Y.%m.%d')
        dir0m = Ldir['LOo'] + 'cast/'
        fnm = dir0m + Ldir['gtagex'] + '/' + station + '_' + date_string + '.nc'
        if testing == True:
            print(fnm)
        try:
            ds = nc4.Dataset(fnm)
            z = ds['z_rho'][:].squeeze()
            z = z - z[-1] # reference to SSH, so top value is always at 0
            for Zn in Z_list:
                i0, i1, fr = zfun.get_interpolant(np.array([Zn]), z, extrap_nan=True)
                for vn in mod_vn_dict.keys():
                    try:
                        mcv = ds[vn][:].squeeze()
                        # note we also apply the scaling factor here
                        val = ((1-fr)*mcv[int(i0)] + fr*mcv[int(i1)])*mod_fac_dict[vn]
                        bc.loc[Zn,mod_vn_dict[vn]] = val
                    except IndexError:
                        pass
            ds.close()
        except OSError:
            if testing == True:
                print('error reading fnm')
            pass
        Bc = pd.concat((Bc,bc), ignore_index=True, sort=False)
            
        
# save the output to disk (only for full set)
if testing == False:
    out_fn = 'ObsMod_' + Ldir['gtagex'] + '_'+ str(year) + '.p'
    print('Saving ' + out_fn)
    Bc.to_pickle(dir11 + out_fn)

# # plotting
fs=16
plt.rc('font', size=fs)
plt.close('all')
# x limits for plotting
dt0 = pd.datetime(year,1,1)
dt1 = pd.datetime(year,12,31)
# y limits for plotting
lim_dict = {'Salinity': (10,35), 'Temp. (deg C)': (5,20),
'DO (mg L-1)': (0, 15), 'DIN (uM)': (0,40)}

# lim_dict = {'Salinity': (0, 34), 'Temp. (deg C)': (0, 24),
#     'DO (mg L-1)': (0, 20), 'DIN (uM)': (0,55), 'Chl (mg m-3)': (0,50)}
lab_dict = {'Salinity':'(a) Salinity', 'Temp. (deg C)':'(b) Temperature [$^{\circ}C$]',
    'DO (mg L-1)': '(c) DO [$mg \ L^{-1}$]', 'DIN (uM)':'(d) DIN  [$\mu M$]'}

if for_web:
    figsize = (6,4)
    NR = 2; NC = 2
else:
    if no_map == True:
        figsize = (15,8)
        NR = 2; NC = 2
    else:
        figsize = (12,5)
        NR = 2; NC = 3
    
clist = ['r', 'g', 'b', 'k']
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
        if (pp==3) and (no_map == False):
            pp = 4
    
    if no_map == False:
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
