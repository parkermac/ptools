#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Look at data by region, binned by month and depth, as monthly
climatologies over selected time periods.

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

for region in range(1,9):
    
    a0 = a[region][0]
    a10 = a[region][1]
    a30 = a[region][2]
    
    count = 0
    am_dict = dict()
    ac_dict = dict()
    y_list = [(1933,1942), (1948,1976), (1999,2017)]
    for y in y_list:
        y0 = pd.datetime(y[0],1,1)
        y1 = pd.datetime(y[1],12,31)
        
        A0 = a0.loc[y0:y1,:]
        A10 = a10.loc[y0:y1,:]
        A30 = a30.loc[y0:y1,:]
        
        ival = np.arange(.5,12.5,1)
        am0 = pd.DataFrame(index=ival, columns=A0.columns)
        am10 = pd.DataFrame(index=ival, columns=A10.columns)
        am30 = pd.DataFrame(index=ival, columns=A30.columns)
        
        # and keep track of count of good points
        ac0 = pd.DataFrame(index=ival, columns=A0.columns)
        ac10 = pd.DataFrame(index=ival, columns=A10.columns)
        ac30 = pd.DataFrame(index=ival, columns=A30.columns)
        
        for m in range(1,13):
            mm = m - 0.5
            
            am0.loc[mm,:] = A0[A0.index.month==m].mean(axis=0)
            am10.loc[mm,:] = A10[A0.index.month==m].mean(axis=0)
            am30.loc[mm,:] = A30[A0.index.month==m].mean(axis=0)

            ac0.loc[mm,:] = A0[A0.index.month==m].count()
            ac10.loc[mm,:] = A10[A0.index.month==m].count()
            ac30.loc[mm,:] = A30[A0.index.month==m].count()
            
        am_dict[count] = (am0,am10,am30)
        ac_dict[count] = (ac0,ac10,ac30)
        count += 1
    
    # plotting
    vn_list = ['Temp. (deg C)', 'Salinity', 'DO (mg L-1)', 'NO3 (uM)']
    lim_list = [[6,24],         [14,34],     [0,14],        [0,34]]
    lim_dict = dict(zip(vn_list, lim_list))
    nn = 1
    fig = plt.figure(figsize=(12,8))
    for vn in vn_list:
        ax = fig.add_subplot(2,2,nn)
        
        alist = [':','--','-']
        for count in range(3):
            am0 = am_dict[count][0]
            am10 = am_dict[count][1]
            am30 = am_dict[count][2]

            ac0 = ac_dict[count][0]
            ac10 = ac_dict[count][1]
            ac30 = ac_dict[count][2]
            
            # drop bins with too few years
            nfew = 3
            am0[ac0<=nfew] = np.nan
            am10[ac10<=nfew] = np.nan
            am30[ac30<=nfew] = np.nan
            
            if (nn==1) and (count==2):
                leg = True
            else:
                leg = False
            am0.plot(y=vn, ax=ax, legend=leg, grid=True, label='Z = 0 m', color='r', linestyle=alist[count])
            am10.plot(y=vn, ax=ax, legend=leg, grid=True, label='Z = -10 m', color='g', linestyle=alist[count])
            am30.plot(y=vn, ax=ax, legend=leg, grid=True, label='Z = -30 m', color='b', linestyle=alist[count])
            ax.set_title(vn)
            ax.set_xlim(0,12)
            ax.set_ylim(lim_dict[vn])
            ax.set_xticks(ival)
            if nn in [3,4]:
                ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
                ax.set_xlabel('Month')
            else:
                ax.set_xticklabels([])
                
            if nn == 4:
                    ax.text(.05,.05,'Solid Lines: Years ' + str(y_list[2][0]) + '-' +str(y_list[2][1]),
                        transform=ax.transAxes)
                    ax.text(.05,.12,'Dashed Lines: Years ' + str(y_list[1][0]) + '-' +str(y_list[1][1]),
                        transform=ax.transAxes)
                    ax.text(.05,.19,'Dotted Lines: Years ' + str(y_list[0][0]) + '-' +str(y_list[0][1]),
                        transform=ax.transAxes)
            
            # if count == 1:
            #     for ii in am0.index:
            #         ax.text(ii,am0.loc[ii,vn],str(ac0.loc[ii,vn]))
            
        nn += 1
    fig.suptitle(collias_region_names[region])
    plt.savefig(dir1 + 'Monthly_' + collias_region_names[region].replace(' ','_') + '.png')

plt.show()