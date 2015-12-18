
"""
Working with EcologyCTD top and bottom time series, and river, and sunshine.
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np

from matplotlib.patches import Rectangle
import matplotlib.dates as mdates

pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart('cascadia1','base')

# Data collections
dc = dict()
dc[0] = ('Hood Canal South', 'HCB003_0', 'HC_HP', (0,8), 'Skokomish', (0,500), (2006,2009)) # GOOD (2006,2014)
dc[1] = ('Hood Canal Twanoh', 'HCB004_0', 'HC_TW', (0,8), 'Skokomish', (0,500), (2005,2009)) # GOOD (2005,2014)
dc[2] = ('Carr Inlet', 'CRR001_0', 'CI', (0,1), 'Nisqually', (0,300), (2011,2013)) # GOOD
dc[3] = ('Main Basin', 'PSB003_0', 'PW', (0,2), 'Green', (0,300), (2010,2011))
dc[4] = ('Hood Canal Dabob', 'HCB010_0', 'HC_DB', (0,6), 'Dosewallips', (0,100), (2010,2014))
dc[5] = ('Hood Canal Duckabush', 'HCB010_0', 'HC_DK', (0,6), 'Duckabush', (0,100), (2006,2009)) # GOOD

##### USER INPUT #######

# Choose the collection
idc = 0

# and choose the time limits
try:
    year_lims = dc[idc][6]
except:
    year_lims = (1998, 2015)
dt0 = datetime(year_lims[0],1,1)
dt1 = datetime(year_lims[1],12,31)

##### END USER INPUT ###

collection = dc[idc][0]
ecy_sta = dc[idc][1]
orca_sta = dc[idc][2]
strat_lims = dc[idc][3]
riv_name = dc[idc][4]
riv_lims = dc[idc][5]

#%% EcologyCTD Data
# These are saved as Series of Sigma (top or bottom layer)

# where are the data csv files
dir0 = '/Users/PM5/Documents/'
indir0 = dir0 + 'ptools_output/env_data/'
indir = indir0 + 'EcologyCTD/'

# choose which stations to include by commenting out parts of this list
# names ending in _0 are multi-year 1990-2014
sta_list_ecy = [
    'BLL009_0', # Bellingham Bay
    'BUD005_0', # Budd Inlet
    'DNA001_0', # Dana Passage
                # 2001 deep T bad in July
                # 2002 shallow S low value Feb, shallow T spike Oct 
                # 2006 lots of little T, s spikes all depths, many months
    'HCB004_0', # Hood Canal in Lynch Cove
    'PSB003_0', # Main Basin off Seattle
                # 1994 salinity spike 4 m Feb/Mar
                # 2007 many T, s spikes, especially deep
    'HCB003_0', # Hood Canal, middle of main channel (Hama Hama)
    'SAR003_0', # middle of Whidbey Basin
    'GOR001_0', # Gordon Point - South Puget Sound
    'PSS019_0', # South Whidbey Basin - off Everett
    'CRR001_0', # Carr Inlet
    'ADM002_0', # Admiralty Inlet just outside Puget Sound
    'ADM003_0', # Admiralty Inlet near Hood Canal
    'HCB010_0', # Hood Canal Near Bangor
    ]

ecy_dft_dict = dict()   
ecy_dfb_dict = dict()      
for sta in sta_list_ecy:
    infile = indir + 'df_top_' + sta + '.p'
    ecy_dft_dict[sta] = pd.read_pickle(infile)
    infile = indir + 'df_bot_' + sta + '.p'
    ecy_dfb_dict[sta] = pd.read_pickle(infile)   
    
#%% ORCA Data
    
# choose station(s) to plot
sta_list_orca = [
    'HC_NB', # North Buoy
    'HC_DB', # Dabob Bay
    'HC_DK', # Duckabush?
    'HC_HP', # Hoodsport
    'HC_TW', # Twanoh
    'PW',    # Pt. Wells
    'CI',    # Carr Inlet
    ]
    
# where are the data csv files
dir0 = '/Users/PM5/Documents/'
indir0 = dir0 + 'ptools_output/env_data/'
indir = indir0 + 'ORCA/'

orca_dft_dict = dict()   
orca_dfb_dict = dict()     
for sta in sta_list_orca:
    infile = indir + 'df_top_' + sta + '.p'
    orca_dft_dict[sta] = pd.read_pickle(infile)
    infile = indir + 'df_bot_' + sta + '.p'
    orca_dfb_dict[sta] = pd.read_pickle(infile)
    
#%% Rivers
       
riv_indir = Ldir['data'] + 'rivers/data_processed/'    
rqt = pd.read_pickle(riv_indir + 'clim_' + riv_name + '.p')

#%% NCEP

dir0 = '/Users/PM5/Documents/'
indir0 = dir0 + 'ptools_output/env_data/'
ncep_indir = indir0 + 'reanalysis/'
ncep_fn = 'df_NCEP_dswrf_Seattle.p'

ncep_df = pd.read_pickle(ncep_indir + ncep_fn)

#%% MM5-WRF
dir0 = '/Users/PM5/Documents/'
indir0 = dir0 + 'ptools_output/env_data/'
mm5_wrf_indir = indir0 + 'mm5_wrf/'
wind_df = pd.read_pickle(mm5_wrf_indir + 'mm5_wrf_jdfE.p')

wind_df['mix'] = (np.abs(wind_df['u10r']**3) + np.abs(wind_df['v10r']**3))

wind_dfm = wind_df.resample('W', how='mean')
for vn in wind_dfm.keys():
    wind_dfm[vn + '_anom'] = wind_dfm[vn]
for vn in wind_df.keys():
    for mo in range(1,13):    
        wind_dfm.ix[wind_dfm.index.month==mo, vn + '_anom'] -= (
            wind_dfm.ix[wind_dfm.index.month==mo, vn].mean(axis=0) )

#%% Plotting set up
    
data_to_plot = [
    #'Salinity',
    #'Temperature',
    #'Sigma',
    'Chl',
    #'DO', 
    ]

# data names, units, and plot ranges
data_names =  ['Salinity','Temperature', 'Sigma',      'Chl',           'Fluor', 'DO']
data_units =  ['psu',     '$^{\circ}C$', '$kg$ $m^{3}$','${\mu}g$ $l^{-1}$','units', '$mg$ $l^{-1}$']
data_ranges = [(14,34),   (4,20),        (14,26),      (0,40),        (0,60),  (0,16)]

# dictionaries to look up data attributes using names
data_unit_dict = dict(zip(data_names, data_units))
data_range_dict = dict(zip(data_names, data_ranges))

# function to add a colored patch to an annual time period
def add_rect(ax, dt0, dt1, y_lims):
    for yr in range(dt0.year, dt1.year + 1):
        startTime = datetime(yr, 2, 1)
        endTime = datetime(yr, 8, 15)
        # convert to matplotlib date representation
        start = mdates.date2num(startTime)
        end = mdates.date2num(endTime)
        width = end - start
        height = y_lims[1] - y_lims[0]
        # Plot rectangle
        rect = Rectangle((start, y_lims[0]), width, height, color='gold', alpha=.3)
        ax.add_patch(rect)

#%% Plotting

plt.close('all')

NR = len(data_to_plot) + 4
NC = 1
fig_size = (20,10)
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=fig_size, squeeze=False)

ecy_dft = ecy_dft_dict[ecy_sta]
ecy_dfb = ecy_dfb_dict[ecy_sta]

orca_dft = orca_dft_dict[orca_sta]
orca_dfb = orca_dfb_dict[orca_sta]

if False:
    orca_dft = orca_dft.resample('W', how='mean', label='left', loffset='3d')
    orca_dfb = orca_dfb.resample('W', how='mean', label='left', loffset='3d')

# Assume the ORCA Fluor is close to Ecology Chl
orca_dft.rename(columns={'Fluor': 'Chl'}, inplace=True)
orca_dfb.rename(columns={'Fluor': 'Chl'}, inplace=True)

ecy_strat = ecy_dfb['Sigma'] - ecy_dft['Sigma']
orca_strat = orca_dfb['Sigma'] - orca_dft['Sigma']

ic = 0
ir = 0

# Ecology and ORCA two-layer data
for fld in data_to_plot:
    
    try:
        axes[ir, ic].plot(ecy_dft.index.values, ecy_dft[fld].values, 'or')
        axes[ir, ic].plot(ecy_dfb.index.values, ecy_dfb[fld].values, 'ob')
    except KeyError:
        pass
    
    try:
        axes[ir, ic].plot(orca_dft.index.values, orca_dft[fld].values, '-r')
        axes[ir, ic].plot(orca_dfb.index.values, orca_dfb[fld].values, '-b')
    except KeyError:
        pass
    
    axes[ir, ic].set_xlim(dt0, dt1)
    axes[ir, ic].set_ylim(data_range_dict[fld])   
    axes[ir, ic].xaxis.set_ticklabels([])
    axes[ir, ic].grid()
    add_rect(axes[ir, ic], dt0, dt1, data_range_dict[fld])
        
    if ir == 0:
        axes[ir, ic].set_title('%s /// EcologyCTD=%s /// ORCA=%s' %
                               (collection, ecy_sta, orca_sta))
        
    if ic == 0:
        axes[ir, ic].text(.05, .8, fld + ' (' + data_unit_dict[fld] + ')', fontsize=16,
                          transform=axes[ir, ic].transAxes)
    ir += 1

# Wind Mixing
axw = axes[ir, ic]
axw.plot(wind_dfm.index.values, wind_dfm['mix'].values)
axw.set_xlim(dt0, dt1)
mix_lims = (0, 2500)
axw.set_ylim(mix_lims)
axw.text(.05, .8, 'Wind Mixing' + ' ($m^{3}$ $s^{-3}$)', fontsize=16,
                  transform=axes[ir, ic].transAxes)
add_rect(axw, dt0, dt1, mix_lims)
ir += 1

# river flow   
axes[ir, ic].plot(rqt.index.values, rqt.values, '-k')
axes[ir, ic].set_xlim(dt0, dt1)
axes[ir, ic].set_ylim(riv_lims)
axes[ir, ic].text(.05, .8, riv_name + ' ($m^{3}$ $s^{-1}$)', fontsize=16,
                  transform=axes[ir, ic].transAxes)
axes[ir, ic].xaxis.set_ticklabels([])
axes[ir, ic].grid()
add_rect(axes[ir, ic], dt0, dt1, riv_lims)
ir += 1

# stratification   
axes[ir, ic].plot(ecy_strat.index.values, ecy_strat.values, 'ok')
axes[ir, ic].plot(orca_strat.index.values, orca_strat.values, '-k')
axes[ir, ic].set_xlim(dt0, dt1)
axes[ir, ic].set_ylim(strat_lims)
axes[ir, ic].text(.05, .8, 'Stratification ($kg$ $m^{-3}$)', fontsize=16,
                  transform=axes[ir, ic].transAxes)
axes[ir, ic].xaxis.set_ticklabels([])
axes[ir, ic].grid()
add_rect(axes[ir, ic], dt0, dt1, strat_lims)
ir += 1
                 
# Atmosphere
x = ncep_df.index.values
y = ncep_df['anom'].values
axes[ir, ic].plot(x, y, '-k')
axes[ir, ic].plot([dt0, dt1], [0, 0], '-k')
d = np.zeros(len(y))
axes[ir, ic].fill_between(x, y, where=y>=d, interpolate=True, color='lightsalmon')
axes[ir, ic].fill_between(x, y, where=y<=d, interpolate=True, color='slateblue')

# add MM5-WRF data
wind_dfm['swdown_anom'].plot(ax=axes[ir, ic])

axes[ir, ic].set_xlim(dt0, dt1)
sw_lims = (-100, 100)
axes[ir, ic].set_ylim(sw_lims)
axes[ir, ic].text(.05, .8, 'NCEP Reanalysis Monthly Downward SW Rad. Anomaly [Seattle] ($W$ $m^{-2}$)', fontsize=16,
                  transform=axes[ir, ic].transAxes)
axes[ir, ic].grid()
add_rect(axes[ir, ic], dt0, dt1, sw_lims)
axes[ir, ic].set_xlabel('Year', fontsize=20)



ir += 1
        
plt.show()