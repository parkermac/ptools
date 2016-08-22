# -*- coding: utf-8 -*-
"""
Looking at anomaly fields for correlation.  LLTK

Created on Sat Dec  5 13:38:08 2015

@author: PM5
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta

pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart('cascadia1','base')

#%% EcologyCTD Data
# These are saved as Series of Sigma (top or bottom layer)

# where are the data csv files
dir0 = '/Users/PM5/Documents/'
indir0 = dir0 + 'ptools_output/env_data/'
indir = indir0 + 'EcologyCTD/'

# choose which stations to include by commenting out parts of this list
# names ending in _0 are multi-year 1990-2014
sta_list_ecy = [
#    'BLL009_0', # Bellingham Bay
#    'BUD005_0', # Budd Inlet
#    'DNA001_0', # Dana Passage
#                # 2001 deep T bad in July
#                # 2002 shallow S low value Feb, shallow T spike Oct 
#                # 2006 lots of little T, s spikes all depths, many months
#    'HCB004_0', # Hood Canal in Lynch Cove
#    'PSB003_0', # Main Basin off Seattle
#                # 1994 salinity spike 4 m Feb/Mar
#                # 2007 many T, s spikes, especially deep
#    'HCB003_0', # Hood Canal, middle of main channel (Hama Hama)
#    'SAR003_0', # middle of Whidbey Basin
#    'GOR001_0', # Gordon Point - South Puget Sound
#    'PSS019_0', # South Whidbey Basin - off Everett
    'CRR001_0', # Carr Inlet
#    'ADM002_0', # Admiralty Inlet just outside Puget Sound
#    'ADM003_0', # Admiralty Inlet near Hood Canal
#    'HCB010_0', # Hood Canal Near Bangor
    ]

ecy_dft_dict = dict()   
ecy_dfb_dict = dict()      
for sta in sta_list_ecy:
    infile = indir + 'df_top_' + sta + '.p'
    ecy_dft_dict[sta] = pd.read_pickle(infile)
    infile = indir + 'df_bot_' + sta + '.p'
    ecy_dfb_dict[sta] = pd.read_pickle(infile)
    
#%% NCEP Sunshine Re-analaysis
    
dir0 = '/Users/PM5/Documents/'
indir0 = dir0 + 'ptools_output/env_data/'
ncep_indir = indir0 + 'reanalysis/'
ncep_fn = 'df_NCEP_dswrf_Seattle.p'

ncep_df = pd.read_pickle(ncep_indir + ncep_fn)

#%% Setup for plotting
    
months = range(1,13) # a generator for a list of 1 to 12

month_color_dict = dict(zip(months,
    ['mediumblue',
    'royalblue', 
    'cadetblue',
    'aquamarine',
    'lightgreen',
    'greenyellow',
    'gold',
    'orange',
    'lightsalmon',
    'mediumorchid',
    'slateblue',
    'purple']))
       
month_name_dict = dict(zip(months,
    ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']))

nrv = [0,0,0,0,1,1,1,1,2,2,2,2]
ncv = [0,1,2,3,0,1,2,3,0,1,2,3]
    
    
#%% Make and plot fields
plt.close('all')

for ecy_sta in sta_list_ecy:
    
    ecy_dft = ecy_dft_dict[ecy_sta]
    ecy_dfb = ecy_dfb_dict[ecy_sta]
    strat = ecy_dfb['Sigma'] - ecy_dft['Sigma'] # a Series
    chl = ecy_dft['Chl']       # a Series
    
    strata = pd.DataFrame(index=strat.index, columns=['fld', 'clim', 'anom'])
    strata['fld'] = strat.values
    strata['clim'] = 0.
    strata['anom'] = strat.values
    for mo in months:
        try:
            mo_mean = strata.ix[strata.index.month==mo, 'fld'].mean(axis=0)
            strata.ix[strata.index.month==mo, 'clim'] = mo_mean
            strata.ix[strata.index.month==mo, 'anom'] -= mo_mean
        except:
            pass
        
    chla = pd.DataFrame(index=chl.index, columns=['fld', 'clim', 'anom'])
    chla['fld'] = chl.values
    chla['clim'] = 0.
    chla['anom'] = chl.values
    for mo in months:
        try:
            mo_mean = chla.ix[chla.index.month==mo, 'fld'].mean(axis=0)
            chla.ix[chla.index.month==mo, 'clim'] = mo_mean
            chla.ix[chla.index.month==mo, 'anom'] -= mo_mean
        except:
            pass

        
    # Plotting       
    NR = 3
    NC = 4
    fig_size = (15,10)
    fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=fig_size, squeeze=False)
    
    to_plot = 'anom'
    
    xx = strata[to_plot]
    yy = chla[to_plot]    
    
    zz = yy.copy()
    counter = 0
    for date in zz.index:
        date1 = date - timedelta(30) # look at sunshine for a month earlier
        zz[counter] = ncep_df.ix[(ncep_df.index.year==date1.year)
                   & (ncep_df.index.month==date1.month), 'anom'].values[0]
        counter += 1
           
    xl = (np.floor(xx.min()), np.ceil(xx.max()))
    yl = (np.floor(yy.min()), np.ceil(yy.max()))
    zl = (np.floor(zz.min()), np.ceil(zz.max()))
    
    plot_type = 'mixed'
    
    for mo in months:
        
        nr = nrv[mo-1]
        nc = ncv[mo-1]
        ax = axes[nr, nc]
        x = strata.ix[strata.index.month==mo, to_plot]
        y = chla.ix[chla.index.month==mo, to_plot]
        z = zz[zz.index.month==mo]
        
        if plot_type == 'mixed':
            # preparing for fitting
            ch = np.array(y.values, dtype=float)
            st = np.array(x.values, dtype=float)
            sw = np.array(z.values, dtype=float)    
            mask = (~np.isnan(st)) & (~np.isnan(ch)) & (~np.isnan(sw))    
            Ch = ch[mask]
            St = st[mask]
            Sw = sw[mask]   
            A = np.vstack([St, Sw, np.ones(len(St))]).T
            m1, m2, c = np.linalg.lstsq(A, Ch)[0]
            xxx = m1*St + m2*Sw + c
            xl = (np.floor(xxx.min()), np.ceil(xxx.max()))
            ax.plot(xxx, Ch, 'o', color=month_color_dict[mo], markersize=15)
            if mo == 1:
                ax.set_title(ecy_sta + ' /// ' + to_plot.upper())
            ax.text(.1, .85, month_name_dict[mo] , transform=ax.transAxes, fontsize=24)
            ax.set_xlim(xl)
            ax.set_ylim(yl)
            if nr == 2:
                ax.set_xlabel('Mixed Model')
            if nc == 0:
                ax.set_ylabel('Chl Anomaly (${\mu}g$ $l^{-1}$)')                
        elif plot_type == 'strat':
            ax.plot(x, y, 'o', color=month_color_dict[mo], markersize=15)
            if mo == 1:
                ax.set_title(ecy_sta + ' /// ' + to_plot.upper())
            ax.text(.1, .85, month_name_dict[mo] , transform=ax.transAxes, fontsize=24)
            ax.set_xlim(xl)
            ax.set_ylim(yl)
            if nr == 2:
                ax.set_xlabel('Stratification Anomaly ($kg$ $m^{3}$)')
            if nc == 0:
                ax.set_ylabel('Chl Anomaly (${\mu}g$ $l^{-1}$)')
        elif plot_type == 'sw':
            ax.plot(z, y, 'o', color=month_color_dict[mo], markersize=15)
            if mo == 1:
                ax.set_title(ecy_sta + ' /// ' + to_plot.upper())
            ax.text(.1, .85, month_name_dict[mo] , transform=ax.transAxes, fontsize=24)
            ax.set_xlim(zl)
            ax.set_ylim(yl)
            if nr == 2:
                ax.set_xlabel('SW Rad Anomaly ($W$ $m^{-2}$)')
            if nc == 0:
                ax.set_ylabel('Chl Anomaly (${\mu}g$ $l^{-1}$)')
        
        ax.grid()

plt.show()

    