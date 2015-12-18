# -*- coding: utf-8 -*-
"""
Plot ORCA processed data
"""

# setup

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# choose station(s) to plot
sta_name_list = [
    'HC_NB',
    'HC_DB',
    'HC_DK',
    'HC_HP',
    'HC_TW',
    'PW',
    'CI',
    ]
    
# choose variables to process
data_to_plot = [
    'Salinity',
    'Temperature',
    'Sigma',
    'Fluor',
    'DO', 
    ]
    
data_names =  ['Salinity', 'Temperature', 'Sigma',      'Fluor', 'DO']
data_units =  ['psu',      '$^{\circ}C$', '$kg m^{3}$', 'units', '$mg l^{-1}$']
data_ranges = [(14,34),    (4,20),        (14,26),      (0,60),  (0,16)]
data_unit_dict = dict(zip(data_names, data_units))
data_range_dict = dict(zip(data_names, data_ranges))
    
# where are the data csv files
dir0 = '/Users/PM5/Documents/'
indir0 = dir0 + 'ptools_output/env_data/'
indir = indir0 + 'ORCA/'

dft_dict = dict()   
dfb_dict = dict()     
for sta in sta_name_list:
    infile = indir + 'df_top_' + sta + '.p'
    dft_dict[sta] = pd.read_pickle(infile)
    infile = indir + 'df_bot_' + sta + '.p'
    dfb_dict[sta] = pd.read_pickle(infile)

#%% Plotting
plt.close('all')

NR = len(sta_name_list)
NC = len(data_to_plot)
fig_size = (25,15)
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=fig_size, squeeze=False)

ir = 0

for sta in sta_name_list:
    dft = dft_dict[sta]
    dfb = dfb_dict[sta]
    
    ic = 0    
    for fld in data_to_plot:
        axes[ir, ic].plot(dft.index.values, dft[fld].values, color='r')
        axes[ir, ic].plot(dfb.index.values, dfb[fld].values, color='b')
        axes[ir, ic].set_xlim(datetime(2005,1,1), datetime(2015,12,31))
        axes[ir, ic].set_ylim(data_range_dict[fld])
        if ir == NR-1:
            # this is a pathetic hack, but it works
            labels = ['','2006','','2008','','2010','','2012','','2014','']
            axes[ir, ic].set_xticklabels(labels)
        else:
            axes[ir, ic].xaxis.set_ticklabels([])
        if ir == 0:
            axes[ir, ic].set_title(fld + ' ' + data_unit_dict[fld])
        if ic == 0:
            axes[ir, ic].text(.05, .8, sta, fontsize=16,
                              transform=axes[ir, ic].transAxes)
        ic += 1
        
    ir += 1

plt.show()