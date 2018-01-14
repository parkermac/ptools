#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 16:21:09 2017

@author: PM5

Code to plot year-long records of many properties
from a single NDBC buoy.

"""

# import modules
import matplotlib.pyplot as plt
from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

from importlib import reload
import ndbc_fun as ndf
reload(ndf)

# **** USER EDITS ****

# choose station(s) (a list of strings) and a year (int)

sn_list = ['46005']; yr = 2003

#sn_list = ['46005', '46041', '46002', '46015', '46059', '46014', '42040', '42035']; yr = 2003

#sn_list = ['sisw1']; yr = 2015
   
# choose time filtering length (a string) e.g.:
#    m = month
#    w = week
#    d = day
#    h = hour (no filtering)
tf = 'd'

# **** END USER EDITS ****

data_dir0 = '../../ptools_data/ndbc/'
name_dict = ndf.get_name_dict()

# process the data

for sn in sn_list:

    id = sn + 'h' + str(yr)
    fn = data_dir0 + sn + '/' + id + '.txt'
    try:        
        dff = ndf.get_data(fn, tf, yr)
    except:
        # missing year
        pass
    
    # PLOTTING
          
    # single plot for each station
    
    fig = plt.figure(figsize = (13, 8))
    
    lw = 2 # line width
    al = .2 # transparency of filled colors
      
    yearday = dff.index.dayofyear
    taux = dff.taux.values
    tauy = dff.tauy.values
    atmp = dff.ATMP.values
    wtmp = dff.WTMP.values
    wvht = dff.WVHT.values
    dpd = dff.DPD.values
    
    ax = fig.add_subplot(221)
    ax.plot(yearday, tauy, '-k',linewidth=1)
    ax.fill_between(yearday, tauy, 0*tauy,
                     where=tauy>0, color='m', interpolate=True, alpha=al)
    ax.fill_between(yearday, tauy, 0*tauy,
                     where=tauy<0, color='g', interpolate=True, alpha=al)
    ax.grid()
    ax.set_xlim(0, 365)
    ax.set_ylim(-.4, .4)
    ax.set_title(name_dict[sn] + ': ' + str(sn) + ' Year=' + str(yr))
    ax.set_ylabel('N-S Windstress [Pa]')
#    ax.set_xticklabels([])
#    ax.set_xlabel('')
    
    ax = fig.add_subplot(223)
    ax.plot(yearday, taux, '-k',linewidth=1)
    ax.fill_between(yearday, taux, 0*taux,
                     where=taux>0, color='m', interpolate=True, alpha=al)
    ax.fill_between(yearday, taux, 0*taux,
                     where=taux<0, color='g', interpolate=True, alpha=al)
    ax.grid()
    ax.set_xlim(0, 365)
    ax.set_ylim(-.4, .4)
    ax.set_xlabel('Yearday')
    ax.set_ylabel('E-W Windstress [Pa]')
    
    ax = fig.add_subplot(222)
    ax.plot(yearday, atmp, '-r', linewidth=lw)
    ax.plot(yearday, wtmp, '-b', linewidth=lw)
    ax.legend(('Air Temperature [degC]', 'Water Temperature [degC]'),
                loc='lower right')
    ax.grid()
    ax.set_xlim(0, 365)
    ax.set_ylim(5, 35)
    #ax.set_xlabel('Yearday')        
    ax.set_title('Filter=' + tf.upper())
#    ax.set_xticklabels([])
#    ax.set_xlabel('')
    
    ax = fig.add_subplot(224)
    ax.plot(yearday, dpd, '-g', linewidth=lw)
    ax.plot(yearday, wvht, '-m', linewidth=lw)
    ax.legend(('Dominant Wave Period [s]', 'Wave Height [m]'),
                loc='upper right')
    ax.grid()
    ax.set_xlim(0, 365)
    ax.set_ylim(0, 20)
    ax.set_xlabel('Yearday')        
           
    plt.show()
    