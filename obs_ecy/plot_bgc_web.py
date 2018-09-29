#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots data from the Department of Ecology, combining data from CTD casts
and bottles at the same station.  Also including model output.

Desigend to create output that can easily be placed on a validation
page on my website.

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

# +++ load ecology CTD cast data +++
dir0 = Ldir['parent'] + 'ptools_data/ecology/'
# load processed station info and data
sta_df = pd.read_pickle(dir0 + 'sta_df.p')
# add Canadian data
dir1 = Ldir['parent'] + 'ptools_data/canada/'
# load processed station info and data
sta_df_ca = pd.read_pickle(dir1 + 'sta_df.p')
sta_df = pd.concat((sta_df, sta_df_ca), sort=True)

# write station info to the screen in a format that can be pasted
# into a webpage javaScript.
print('sta_list = [')
for station in sta_df.index:
    lon = sta_df.loc[station, 'Longitude']
    lat = sta_df.loc[station, 'Latitude']
    long_name = sta_df.loc[station, 'Descrip']
    if 'SOG' in station:
        long_name = 'Strait of Georgia'
    print('    {sta:"%s", lat:%0.4f, lng:%0.4f, longname:"%s"},' % (station, lat, lon, long_name))
print('];')

year = 2017
Casts = pd.read_pickle(dir0 + 'Casts_' + str(year) + '.p')
Casts_ca = pd.read_pickle(dir1 + 'Casts_' + str(year) + '.p')
Casts = pd.concat((Casts, Casts_ca))

# +++ load ecology bottle data +++
Bottles = pd.read_pickle(dir0 + 'Bottles_' + str(year) + '.p')

Ldir['gtagex'] = 'cas4_v2_lo6biom'

testing = False

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
    dir1 = dir11 + 'bgc_series_web' + str(year) + '_' + Ldir['gtagex'] + '/'
    Lfun.make_dir(dir1, clean=True)


Bc = pd.DataFrame()
for station in sta_to_plot:
    
    # loop over all casts at this station
    #print(' - Working on: ' + station)
    
    casts = Casts[Casts['Station'] == station]   
    casts = casts.set_index('Date')    
    # identify a single cast by its DATE
    calldates = casts.index
    castdates = calldates.unique() # a short list of unique dates (1 per cast)
    
    bottles = Bottles[Bottles['Station'] == station]   
    bottles = bottles.set_index('Date')    
    # identify a single cast by its DATE
    balldates = bottles.index
    bottledates = balldates.unique() # a short list of unique dates (1 per cast)
    
    dates = bottledates.union(castdates)
    
    # some useful dictionaries for renameing and rescaling
    cast_vn_dict = {'Salinity': 'Salinity', 'Temperature': 'Temp. (deg C)',
        'Sigma': 'Sigma (kg m-3)', 'Chl': 'Chlorophyl', 'DO': 'DO (mg L-1)',
        'Turb': 'Turbidity'}
        
    bot_vn_dict = {'PO4(uM)D': 'PO4 (uM)', 'SiOH4(uM)D': 'SiOH4 (uM)',
        'NO3(uM)D': 'NO3 (uM)', 'NO2(uM)D': 'NO2 (uM)', 'NH4(uM)D': 'NH4 (uM)',
       'DIN': 'DIN (uM)'}
       
    mod_fac_dict = {'salt': 1, 'temp': 1,
        'NO3': 1, 'phytoplankton': 1, 'zooplankton': 1,
        'detritus': 1, 'Ldetritus': 1, 'oxygen': 0.032,
        'TIC': 1, 'alkalinity': 1}
    
    # these have model units AFTER multiplication by mod_fac_dict
    mod_vn_dict = {'salt': 'Mod Salinity', 'temp': 'Mod Temp. (deg C)',
        'NO3': 'Mod DIN (uM)', 'phytoplankton': 'Mod Chl', 'zooplankton': 'Mod Zoop',
        'detritus': 'Mod Detr.', 'Ldetritus': 'Mod Lg. Detr.', 'oxygen': 'Mod DO (mg L-1)',
        'TIC': 'Mod DIC (uM C)', 'alkalinity': 'Mod Alkalinity (uEq. L-1)'}
        
    casts = casts.rename(columns=cast_vn_dict)
    bottles = bottles.rename(columns=bot_vn_dict)
        
    for dd in dates:
        # NOTE the brackets around [dd] keep the result as a DataFrame even if
        # we are only pulling out a single row.
        ca = casts.loc[[dd],:]
        try:
            bo = bottles.loc[[dd],:]
        except KeyError:
            bo  = pd.DataFrame.from_dict({'Z':[0,-10,-30], 'Station':station, 'Znom':[0,-10,-30]})
        ca = ca.set_index('Z')
        bo = bo.set_index('Z')
        bc = bo.copy() # a DataFrame to store both bottles and casts, but only at shared depths
        
        for vn in cast_vn_dict.values():
            bc[vn] = np.nan
        
        for Z in bo.index:
            try:
                i0, i1, fr = zfun.get_interpolant(np.array([Z]), ca.index.values, extrap_nan=False)
                for vn in cast_vn_dict.values():
                    cv = ca[vn]
                    bc.loc[Z,vn] = (1-fr)*cv.iloc[int(i0)] + fr*cv.iloc[int(i1)]
            except:
                pass
        bc['Date'] = dd
        bc['Z'] = bc.index
        
        # also get the model fields for this day
        date_string = dd.strftime('%Y.%m.%d')
        dir0m = Ldir['LOo'] + 'casts/'
        fnm = dir0m + Ldir['gtagex'] + '/' + station + '_' + date_string + '.nc'
        try:
            ds = nc4.Dataset(fnm)
            z = ds['z_rho'][:].squeeze()
            z = z - z[-1] # reference to SSH
            
            for Z in bo.index:
                i0, i1, fr = zfun.get_interpolant(np.array([Z]), z, extrap_nan=False)
                for vn in mod_vn_dict.keys():
                    mcv = ds[vn][:].squeeze()
                    val = ((1-fr)*mcv[int(i0)] + fr*mcv[int(i1)])*mod_fac_dict[vn]
                    bc.loc[Z,mod_vn_dict[vn]] = val
            ds.close()
        except OSError:
            pass
        Bc = pd.concat((Bc,bc),ignore_index=True, sort=False)
        
# some statistics
print('')
for vn in ['Temp. (deg C)', 'Salinity', 'DO (mg L-1)', 'DIN (uM)']:
    rmse = np.sqrt( ((Bc[vn]-Bc['Mod '+vn])**2).mean() )
    print('RMSE %s = %0.2f' % (vn, rmse))

# # plotting
plt.close('all')
# x limits for plotting
dt0 = pd.datetime(year,1,1)
dt1 = pd.datetime(year,12,31)
# y limits for plotting
lim_dict = {'Salinity': (0, 34), 'Temp. (deg C)': (0, 24),
    #'DO (mg L-1)': (0, 14), 'DIN (uM)': (0,35)}
    'DO (mg L-1)': (0, 20), 'DIN (uM)': (0,55)}
#print('\nMAKING PLOTS\n')

for station in sta_to_plot:
    fig = plt.figure(figsize=(6,4))
    
    # loop over all casts at this station
    #print(' - Plotting: ' + station)
    
    A = Bc[Bc['Station'] == station]   
    A = A.set_index('Date')
    clist = ['r', 'g', 'b']
    pp = 1
    for vn in ['Salinity', 'Temp. (deg C)', 'DO (mg L-1)', 'DIN (uM)']:
        ax = fig.add_subplot(2,2,pp)
        ii = 0
        for zz in [0, -10, -30]:
            a = A[A['Znom'] == zz]
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
        
        ax.text(.05, .9, vn, fontweight='bold', transform=ax.transAxes)
        if pp==1:
            ax.text(.05, .4, 'Z =   0 m', color='r', fontweight='bold', transform=ax.transAxes)
            ax.text(.05, .3, 'Z = -10 m', color='g', fontweight='bold', transform=ax.transAxes)
            ax.text(.05, .2, 'Z = -30 m', color='b', fontweight='bold', transform=ax.transAxes)
            ax.text(.05, .1, 'Solid = Observed, Dashed = Model',
                fontweight='bold', transform=ax.transAxes)
        if pp in [1,2]:
            ax.set_xticklabels('')
            ax.set_xlabel('')
        else:
            ax.set_xlabel('Date ' + str(year))
        ax.grid(True)
        pp += 1
            
        
    fig.tight_layout()
    
    if save_fig:
        plt.savefig(dir1 + station + '.png')
        plt.close()
    else:
        plt.show()
