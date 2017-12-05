#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 14:16:48 2017

@author: PM5

Plot data on river flow and nearby NDBC buoy wave data.

"""

# import modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

from importlib import reload
import ndbc_fun as ndf
reload(ndf)

data_dir0 = '../../ptools_data/ndbc/'
name_dict = ndf.get_name_dict()

do_save = True
if do_save:
    out_dir00 = '../../ptools_output/'
    out_dir0 = out_dir00 + 'ndbc/'
    ndf.make_dir(out_dir00)
    ndf.make_dir(out_dir0)

# **** USER EDITS ****

# USER uncomment a single sn_list (remove the #)

# the original 3 pairs from the west coast
#sn_list = ['46005', '46041', '46002', '46015', '46059', '46014']

# a single station
#sn_list = ['42040']

# two stations from the Gulf Coast
#sn_list = ['42040', '42035']

# several
#sn_list = ['46005', '46041', '46002', '46015', '46059', '46014', '42040', '42035','sisw1', 'wpow1']

# one station
sn_list = ['46002']

# **** END USER EDITS ****

name_dict = ndf.get_name_dict()

# set time limits to plot

if True: # plot the full record of data I downloaded
    y0 = 1970 # start at the beginning of this year
    y1 = 2020 # go the the end of this year
else: # or select a smaller time span, like a year
    y0 = 2015
    y1 = 2016
    
t0 = datetime(y0,1,1)
t1 = datetime(y1,12,31)

# always choose monthly time filtering
tf = 'm'

for sn in sn_list:
    
    out_fn = out_dir0 + sn + '_monthly.csv'
    out_fnp = out_dir0 + sn + '_monthly.p'

    # # initialize the pandas DataFrame
    # nyears = y1 - y0 + 1
    # iii = ( pd.date_range(str(y0) + '-1-1', periods=12*nyears, freq='M')
    #     - timedelta(days=15) )
    # DFF = pd.DataFrame(index=iii)

    # process all the NDBC data
    yr_list = range(y0,y1+1)
    id_list = []
    first_good_year = True
    for yr in yr_list:
        id = str(sn) + 'h' + str(yr)
        id_list.append(id)
        fn = (data_dir0 + sn + '/' + id + '.txt')
        try:
            dff = ndf.get_data(fn, tf, yr)
            if first_good_year:
                DFF = dff.copy()
                first_good_year = False
            else:
                DFF = DFF.combine_first(dff)
            print(fn + ' = success')
            
            # create a dict relating variable names to their units
            header = pd.read_csv(fn, nrows=1, delim_whitespace=True)
            unit_dict = dict(header.iloc[0])
            unit_dict['taux'] = 'Pa'
            unit_dict['tauy'] = 'Pa'
            unit_dict['wspd_clim'] = unit_dict['WSPD']
            # we re-create it each time because only later years have
            # the units in the form we are looking for
            
        except:
            print(fn + ' = fail')
            pass

    #%% make climatology
    if False:
        DFF['wspd_clim'] = np.nan
        for t in DFF.index:
            mo = t.month
            the_mean = DFF.loc[DFF.index.month==mo, 'WSPD'].mean(axis=0)
            DFF.loc[t, 'wspd_clim'] = the_mean

    #%% PLOTTING

    fig1 = plt.figure(figsize = (14, 8))
    NR = 2
    ax1 = fig1.add_subplot(NR,1,1)
    ax2 = fig1.add_subplot(NR,1,2)
    #ax3 = fig1.add_subplot(NR,1,3)

    ax = ax1
    vn = 'tauy'
    DFF[vn].plot(ax=ax, grid=True, xlim=(t0,t1), ylim=(-.2, .2), color='r')
    ax.set_xticklabels([])
    ax.set_xlabel('')
    ax.set_ylabel(vn + ' [' + unit_dict[vn] + ']')
    ax.set_title('Station ' + str(sn) + ' ' + name_dict[sn])

    ax = ax2
    vn1 = 'ATMP'
    DFF[vn1].plot(ax=ax, grid=True, xlim=(t0,t1), ylim=(0, 20),color='r')
    vn2 = 'WTMP'
    DFF[vn2].plot(ax=ax, grid=True, xlim=(t0,t1), ylim=(0, 20),color='b')
    ax.legend(('Air Temperature [degC]', 'Water Temperature [degC]'),
                loc='upper left')    
    #ax.set_ylabel(vn2 + ' [' + unit_dict[vn2] + ']')

    plt.show()
    
    if do_save:
        DFF.to_csv(out_fn)
        DFF.to_pickle(out_fnp)
    
