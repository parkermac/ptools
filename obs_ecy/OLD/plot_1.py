"""
Plots data from a Department of Ecology csv file, for flight CTD stations.

"""

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np


dir0 = '/Users/pm7/Documents/ptools_data/ecology/'

#sta = 'HCB010_0_provisional'
sta = 'HCB008_0'

# a function to read a csv file into a "data frame"
def get_casts(sta, dir0):
    filename = sta + '.csv'
    filename2 = sta + '_fixed.csv'
    fn = dir0 + filename
    fn2 = dir0 + filename2
    # replace bad data values
    # (hard to parse a csv file is some fields have commas!)
    ff = open(fn, 'r')
    ff2 = open(fn2, 'w')
    rep_list = ['-9,999.0', '(9999.0)', '(9999.00)', '9999.0']
    for line in ff:
        for rep in rep_list:
            line = line.replace(rep,'9999.0')
        ff2.write(line)
    ff.close()
    ff2.close()
    # read csv data into a data frame
    casts = pd.read_csv(fn2, parse_dates = [' date'])
    # remove spaces from column headings
    cols = casts.columns
    cols2 = []
    for col in cols:
        Col = col.replace(' ','')
        cols2.append(Col)
    casts.columns = cols2
    # and specify the index column
    casts = casts.set_index('date')
    # and make a Z column
    casts['Z'] = -casts['depth(meters)']
    # and replace missing data with NaN
    casts[casts==9999] = np.nan
    return casts
    
months = range(1,13) # a list of 1 to 12

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

    
casts = get_casts(sta, dir0)

# Index(['station', 'time', 'depth(meters)', 'salinity(psu)',
#        'temperature(centigrade)', 'density(sigmat)', 'chlorophyllraw(ug/l)',
#        'dissolvedoxygen(mg/l)raw', 'dissolvedoxygen(mg/l)adjusted',
#        'lighttransmission(%)', 'pH', 'Z'],

# identify a single cast by its date
alldates = casts.index
castdates0 = alldates.unique() # a short list of unique dates (1 per cast)
# get a list of all years with casts
years = castdates0.year.unique()

    
# plotting
plt.close('all')

#vn_list = ['salinity(psu)', 'temperature(centigrade)', 'lighttransmission(%)']
vn_list = ['density(sigmat)']
NC = len(vn_list)
NR = 1

for fyear in years:
    
    # select a single year
    dt0 = datetime(fyear,1,1)
    dt1 = datetime(fyear,12,31)
    castdates = castdates0[castdates0 >= dt0]
    castdates = castdates[castdates <= dt1]
    
    fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(12,8), squeeze=False)

    # plot the CTD cast data for this station
    flag = 1
    for cdate in castdates:
        imo = cdate.month
        cast = casts[casts.index==cdate]
        # drop repeat values (aiming for the first of a depth pair)
        zdf = np.diff(cast['Z'])
        zdf = np.concatenate((np.array([1.,]),zdf))
        mask = zdf != 0
        cast = cast[mask]
        cast = cast[:-5] # drop bottom values (sometimes bad)
    
        ncol = 0
        for vn in vn_list:
            ax = axes[0,ncol]
    
            ax.plot(cast[vn].values, cast['Z'],'-', linewidth = 3,
                color=month_color_dict[imo])
            
            if ncol==0 and flag == 1:
                for imo in months:
                    ax.text(.05, 1 - imo/13.,
                        month_name_dict[imo], color=month_color_dict[imo],
                        fontsize=14, fontweight='bold',
                        verticalalignment='center', transform=ax.transAxes)
            ncol += 1
        
            if flag == 1:
                ax.grid()
                ax.set_xlabel(vn)
        flag = 0
    
    axes[0,0].set_title(sta + ': ' + str(fyear))

plt.show()
