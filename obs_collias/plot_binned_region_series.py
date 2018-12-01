#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Look at data by region, binned by month and depth, as long time series.

"""

# imports
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt

collias_region_names = {1:'Strait of Juan de Fuca',
        2:'Admiralty Inlet',
        3:'Main Basin',
        4:'South Sound',
        5:'Hood Canal',
        6:'Whidbey Basin',
        7:'North Sound',
        8:'San Juan Island Passages'}

# directory for processed data to use
dir1 = '../../ptools_output/collias/'

a = pickle.load(open(dir1 + 'region_month_z_means.p', 'rb'))

# plotting
plt.close('all')
vn_list = ['Temp. (deg C)', 'Salinity', 'DO (mg L-1)', 'NO3 (uM)']
lim_list = [[6,24],         [14,34],     [0,14],        [0,34]]
lim_dict = dict(zip(vn_list, lim_list))
for region in range(1,9):
    fig = plt.figure(figsize=(12,12))
    a0 = a[region][0]
    a10 = a[region][1]
    a30 = a[region][2]
    N = len(vn_list)
    ii = 1
    for vn in vn_list:
        ax = fig.add_subplot(N,1,ii)
        if ii == N:
            leg = True
        else:
            leg = False
        a0.plot(y=vn, ax=ax, grid=True, label='Z = 0 m', color='r', legend=leg)
        a10.plot(y=vn, ax=ax, grid=True, label='Z = -10 m', color='g', legend=leg)
        a30.plot(y=vn, ax=ax, grid=True, label='Z = -30 m', color='b', legend=leg)
        ax.set_xlim(pd.datetime(1932,1,1),pd.datetime(2018,1,1))
        if vn=='Salinity':
            ax.text(.05,.15,str(vn),transform=ax.transAxes, fontweight='bold')
        else:
            ax.text(.05,.85,str(vn),transform=ax.transAxes, fontweight='bold')
        ax.set_ylim(lim_dict[vn])
        if ii == 1:
            ax.set_title('Combined Collias and Ecology Data: ' + collias_region_names[region])
        if ii<N:
            ax.set_xticklabels([])
            ax.set_xlabel('')
        if ii == N:
            ax.legend(loc=(.2,.1))
        ii+=1
        
    plt.savefig(dir1 + 'Series_' + collias_region_names[region].replace(' ','_') + '.png')

plt.show()
    
