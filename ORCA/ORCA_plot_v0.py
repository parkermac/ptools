"""
Code to plot ORCA data files.

Parker MacCready
"""

##### USER INPUT ########################

# choose station(s) to plot
sta_name_list = [
    'HC_NB',
    'HC_DB',
    'HC_DK',
    'HC_HP',
    'HC_TW',
    'PW',
    'CI'
    ]

# choose variables to plot
data_to_plot = [
    'Salt',
    'Temp',
    #'Sigma',
    'Fluor',
    'DO',
    #'NO3'  # unreliable, don't use  
    ]

# choose type of plot or output    
do_time_series = True
do_property_property = False
do_csv = False
z_for_csv = -50

# choose z levels to plot    
z_to_plot = [-3, -10, -20, -50, -100]
z_colors = ['r','orange','g','b','k']

# choose start and end times to plot (years - can be decimals)
limit_year = False
year_lims = [2008, 2009]

# station names available : and their long names
sta_name_dict = {'CI':'South Sound - Carr Inlet', 'HC_DB':'Hood Canal - Dabob Bay',
    'HC_DK':'Hood Canal - Duckabush', 'HC_HP':'Hood Canal - Hoodsport',
    'HC_NB':'Hood Canal - North Buoy', 'HC_TW':'Hood Canal - Twanoh',
    'PW':'Main Basin - Point Wells'}

# data location (the path to wherever you have the data files)
# NOTE: Windows users use double forward slashes
# e.g dir0 = 'C://Users//Julie//Downloads//'
dir0 = '/Users/PM5/Documents/tools_data/obs_data/mooring/ORCA_new_2015_02/'

# output location
out_dir = '/Users/PM5/Documents/tools_output/pydev_out/ORCA_out/'

# extra module location (path to where you keep matfun
path_to_alpha = '/Users/PM5/Documents/LiveOcean/alpha'

# data names, units, and plot ranges
data_long_names=['Bsal','Btemp',       'Bsigma',     'Bfluor','Boxy_mgL',    'Bnitrate']
data_names =  ['Salt',  'Temp',        'Sigma',      'Fluor', 'DO',          'NO3']
data_units =  ['psu',   '$^{\circ}C$', '$kg m^{3}$', 'units', '$mg l^{-1}$', '$\mu M$']
data_ranges = [(20,31), (7,20),        (14,26),      (0,80),  (1,14),        (0,40)]

##### END USER INPUT ####################

# setup
import os
import sys
alp = os.path.abspath(path_to_alpha)
if alp not in sys.path:
    sys.path.append(alp)
from importlib import reload
import matfun
reload(matfun)

import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

plt.close('all')

# suppress warning messages (e.g. from YY>Ylo where YY is an array with nan's)
import warnings
warnings.simplefilter("ignore")

# dictionaries to look up data attributes using names
data_name_dict = dict(zip(data_long_names, data_names))
data_unit_dict = dict(zip(data_names, data_units))
data_range_dict = dict(zip(data_names, data_ranges))
z_color_dict = dict(zip(z_to_plot, z_colors))

for sta_name in sta_name_list:
    
    # load the data
    filename = sta_name + '_CTD_data_bin_web'
    fn = dir0 + filename + '.mat'
    a = matfun.loadmat(fn)
    A = dict()
    Z = a['Bpress'] # packed top to bottom
    ZZ = np.nanmean(Z, axis = 1)
    ZZind = np.arange(len(ZZ))
    T = a['Btime']
    NZ, NT = Z.shape
    for long_name in data_name_dict:
        var_name = data_name_dict[long_name]
        A[var_name] = a[long_name]
        A[var_name][A[var_name]==-555.] = np.nan
    
    # make a time axis
    t0 = datetime(2000,1,1)
    T0 = np.nanmin(T)
    delt = timedelta(T0)
    t00 = t0 + delt # datetime of the first good record
    this_year = t00.year # year of the first good record
    t1 = datetime(this_year,1,1)
    ndays = t1.toordinal() - t0.toordinal()
    T = T - ndays
    Y = T/365.25 + this_year
    YY = np.nanmean(Y, axis = 0)
    YYind = np.arange(len(YY))
    
    if limit_year == False: # set time limits automatically
        Ylo = np.floor(np.nanmin(Y))
        Yhi = np.ceil(np.nanmax(Y))
    else: # set time limits by hand
        Ylo = year_lims[0]
        Yhi = year_lims[1]
        
    try:
        iylo = YYind[YY>Ylo].tolist()[0]
        iyhi = YYind[YY<Yhi].tolist()[-1]
    except:
        print('\n' + sta_name)
        print('** This station may not have data for the requested time span. **')
        print('** Plotting full time span instead. **')
        Ylo = np.floor(np.nanmin(Y))
        Yhi = np.ceil(np.nanmax(Y))
        iylo = YYind[YY>Ylo].tolist()[0]
        iyhi = YYind[YY<Yhi].tolist()[-1]
    
    if do_time_series:
        # TIME SERIES PLOTS
        
        # set up the plot axes (an array)
        NR = len(data_to_plot)
        NC = 1
        fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(15,15), squeeze=False)
        
        iir = 0
        for dtp in data_to_plot:
            for zz in z_to_plot:
                try:
                    iz = ZZind[ZZ==zz].tolist()[0]
                    axes[iir, 0].plot(Y[iz,iylo:iyhi], A[dtp][iz,iylo:iyhi], '-', color=z_color_dict[zz]) 
                except:
                    pass  
            axes[iir, 0].set_ylim(data_range_dict[dtp])
            axes[iir, 0].set_xlim((Ylo, Yhi)) # time range to plot
            axes[iir, 0].ticklabel_format(useOffset=False, axis='x')
            axes[iir, 0].set_ylabel(dtp + ' ' + data_unit_dict[dtp])
            for yy in range(int(Ylo) + 1, int(Yhi)):
                axes[iir, 0].plot([yy, yy], data_range_dict[dtp], '-k')
            if iir == 0:
                axes[iir, 0].set_title(sta_name_dict[sta_name])
                zcount = 0
                for zz in z_to_plot:
                    try:
                        axes[iir, 0].text(.01, .4 - .09*zcount, str(zz) + ' m', color=z_color_dict[zz],
                            transform=axes[iir, 0].transAxes, fontsize = 12)
                        zcount += 1
                    except:
                        pass
            if iir == NR-1:
                axes[iir, 0].set_xlabel('YEAR')
            iir += 1
            
        fig.show()
    
    if do_property_property:    
        # PROPERTY-PROPERTY PLOTS       
        # set up the plot axes (an array)
        NR = 1
        NC = 2
        fig2, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(15,8), squeeze=False)
        
        iic = 0
        pairs_to_plot = [('Salt', 'Temp'), ('Fluor', 'DO')]
        for pair in pairs_to_plot:
            for zz in z_to_plot:
                try:
                    iz = ZZind[ZZ==zz].tolist()[0]
                    axes[0, iic].plot(A[pair[0]][iz,iylo:iyhi],
                        A[pair[1]][iz,iylo:iyhi], '.', color=z_color_dict[zz]) 
                except:
                    pass
            axes[0, iic].set_xlim(data_range_dict[pair[0]])
            axes[0, iic].set_ylim(data_range_dict[pair[1]])
            axes[0, iic].set_xlabel(pair[0] + ' ' + data_unit_dict[pair[0]])
            axes[0, iic].set_ylabel(pair[1] + ' ' + data_unit_dict[pair[1]])
            if iic == 0:
                axes[0, iic].set_title(sta_name_dict[sta_name])
                zcount = 0
                for zz in z_to_plot:
                    try:
                        axes[0, iic].text(.01, .9 - .1*zcount, str(zz) + ' m',
                            color=z_color_dict[zz],
                            transform=axes[0, iic].transAxes, fontsize = 18)
                        zcount += 1
                    except:
                        pass
            iic += 1
        
        fig2.show()
        
    if do_csv:
        # output one selected depth level to a csv file
        try:
            iz = ZZind[ZZ==z_for_csv].tolist()[0]
            csv_name = out_dir + 'ORCA_' + sta_name + '_' + str(abs(z_for_csv)) + '.csv'
            write_mode='wb'
            import csv
            with open(csv_name, write_mode) as ff:
                ww = csv.writer(ff)
                ww.writerow(['Date'] + data_to_plot)
                for ii in range(iylo,iyhi):
                    this_row = ['%.4f' % YY[ii]]
                    for dtp in data_to_plot:                  
                        this_row.append('%.4f' % A[dtp][iz,ii])
                    if 'nan' not in this_row: # omit rows with ANY bad data
                        ww.writerow(this_row)
        except:
            print('\n' + sta_name + ' request to make csv output failed')
            print('** This station may not have data at the requested depth **')
        
