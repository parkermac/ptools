# -*- coding: utf-8 -*-
"""
Code to plot year-long records of atmosphere and ocean surface properties
from NDBC buoys.

"""

# import modules
import matplotlib.pyplot as plt
from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages
from importlib import reload
import coast_fun
reload(coast_fun)

# **** USER EDITS ****

# the code will loop over all the regions in region_list

region_list = ['WA', 'OR', 'CA']

for region in region_list:

    # specify stations and year
    if region == 'WA':
        sn_list = ['46005', '46041'] # Washington
        yr = 2003
    elif region == 'OR':
        sn_list = ['46002', '46015'] # Oregon
        yr = 2007
    elif region == 'CA':
        sn_list = ['46059', '46014'] # California
        yr = 2007
    
    # choose time filtering length (a string) e.g.:
    #    m = month
    #    w = week
    #    d = day
    #    h = hour (no filtering)
    tf = 'w' # 'w' looks good
        
    # **** END USER EDITS ****

    name_dict = coast_fun.get_name_dict()
    
    # initialize a dict to organize the output
    dff_dict = dict()
    # process the data
    for sn in sn_list:
        id = sn + 'h' + str(yr)
        fn = sn + '/' + id + '.txt'
        try:        
            dff = coast_fun.get_data(fn, tf, yr)
            dff_dict[id] = dff
        except:
            # missing year
            pass
    
    # PLOTTING
    
    # single plot to compare stations
    
    fig = plt.figure(figsize = (14, 7))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
     
    counter = 0
    for sn in sn_list:
        
        id = str(sn) + 'h' + str(yr)
        
        dff = dff_dict[id]
           
        yearday = dff.index.dayofyear    
        tauy = dff.tauy.values
        atmp = dff.ATMP.values
        wtmp = dff.WTMP.values
        
        lw = 2 # line width
        al = .2 # transparency of filled colors
        
        if counter==0:
            ax1.plot(yearday, tauy, '-k',linewidth=lw)
            ax1.fill_between(yearday, tauy, 0*tauy,
                             where=tauy>0, color='m', interpolate=True, alpha=al)
            ax1.fill_between(yearday, tauy, 0*tauy,
                             where=tauy<0, color='g', interpolate=True, alpha=al)
            ax1.grid()
            ax1.set_xlim(0, 365)
            ax1.set_ylim(-.2, .2)
            ax1.set_title('Year ' + str(yr))
            ax1.set_xlabel('Yearday')
            ax1.text(10,-.15,'UPWELLING', color='g', fontweight='bold')
            ax1.text(10,.15,'DOWNWELLING', color='m', fontweight='bold')
            ax1.set_ylabel('N-S Windstress [Pa]')
        else:
            ax1.plot(yearday, tauy, '--k',linewidth=lw)
            ax1.fill_between(yearday, tauy, 0*tauy,
                             where=tauy>0, color='m', interpolate=True, alpha=al)
            ax1.fill_between(yearday, tauy, 0*tauy,
                             where=tauy<0, color='g', interpolate=True, alpha=al)
    
        if counter==0:
            ax2.plot(yearday, atmp, '-r', linewidth=lw)
            ax2.plot(yearday, wtmp, '-b', linewidth=lw)
            ax2.legend(('Air Temperature [degC]', 'Water Temperature [degC]'),
                        loc='lower right')
            ax2.grid()
            ax2.set_xlim(0, 365)
            ax2.set_ylim(5, 20)
            ax2.set_xlabel('Yearday')        
            ax2.text(.05, .95, 'Solid: ' + name_dict[sn] + ' ' + str(sn),
                transform=ax2.transAxes)
            
        else:
            ax2.plot(yearday, atmp, '--r', linewidth=lw)
            ax2.plot(yearday, wtmp, '--b', linewidth=lw)
            ax2.text(.05, .90, 'Dashed: ' + name_dict[sn] + ' ' + str(sn),
                transform=ax2.transAxes)
       
        counter += 1
        
    plt.show()
    