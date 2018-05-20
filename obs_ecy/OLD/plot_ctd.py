"""
Plots data from a Department of Ecology csv file, for flight CTD stations.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

## USER INPUT ##

# where are the data csv files
dir0 = '/Users/pm7/Documents/ptools_data/ecology/'

# choose which station to plot
# names ending in _0 are multi-year 1990-
sta = 'HCB010_0_provisional'

# choose which data fields to plot by commenting out parts of this list
data_to_plot = [
    'Salinity',
    'Temperature',
    #'Sigma',
    #'Chl',
    #'DO',
    #'Trans',
    #'pH'
    ]

year = 2016 # integer year, or None to plot all years

## END USER INPUT ##

# lists of data properties

# data long names (used in the csv after removing spaces)
# and we retain only these fields
data_long_names = ['salinity(psu)', 'temperature(centigrade)',
                   'density(sigmat)', 'chlorophyllraw(ug/l)',
                   'dissolvedoxygen(mg/l)adjusted',
                   'lighttransmission(%)', 'pH', 'Z']

# data short names, units, and plot ranges
data_names =  ['Salinity','Temperature','Sigma', 'Chl', 'DO',   'Trans',  'pH', 'Z']
data_units =  ['psu',     'deg C',      'kg/m3', 'ug/l', 'mg/l', '%',     '',   'm']
data_ranges = [(14,34),   (4,20),       (14,26), (0,40), (0,18), (0,100), (6,9),(-50,0)]

# dictionaries to look up data attributes using names
data_name_dict = dict(zip(data_long_names, data_names))
data_unit_dict = dict(zip(data_names, data_units))
data_range_dict = dict(zip(data_names, data_ranges))

# more useful lists and dictionaries for plotting different months

months = range(1,13) # a list of 1 to 12

month_color_dict = dict(zip(months,
    ['mediumblue', 'royalblue', 'cadetblue', 'aquamarine',
    'lightgreen', 'greenyellow', 'gold', 'orange',
    'lightsalmon', 'mediumorchid', 'slateblue', 'purple']))

month_name_dict = dict(zip(months,
    ['Jan','Feb','Mar','Apr','May','Jun',
     'Jul','Aug','Sep','Oct','Nov','Dec']))

def make_fixed(sta, dir0):
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
    return fn2

# a function to read a csv file into a "data frame"
def get_casts(fn):
    # read csv data into a data frame
    casts = pd.read_csv(fn, parse_dates = [' date'])
    # remove spaces from column headings
    cols = casts.columns
    cols2 = []
    for col in cols:
        cols2.append(col.replace(' ',''))
    casts.columns = cols2
    casts = casts.set_index('date') # and specify the index column
    casts['Z'] = -casts['depth(meters)'] # and make a Z column
    # and replace missing data with NaN
    casts[casts==9999] = np.nan
    return casts

# PLOTTING

# setup
plt.close('all')
figsize = (13,7)

# plot CASTS

# set up the plot axes (an array)
NR = 1
NC = len(data_to_plot)

fig, axes = plt.subplots(nrows=NR, ncols=NC, sharey=True,
                         figsize=figsize, squeeze=False)

# dictionaries linking data field to column number
get_nc = dict(zip(data_to_plot, range(NC)))


# get this station's data frame using our function
fn = make_fixed(sta, dir0)    
casts = get_casts(fn)
casts = casts[data_long_names] # keep only selected columns
casts = casts.rename(columns=data_name_dict) # and rename them

# identify a single cast by its DATE
alldates = casts.index
castdates = alldates.unique() # a short list of unique dates (1 per cast)

title = sta
if year != None:
    title = sta + ': ' + str(year)
    dt0 = datetime(year,1,1)
    dt1 = datetime(year,12,31)
    castdates = castdates[castdates >= dt0]
    castdates = castdates[castdates <= dt1] 
    casts = casts.loc[dt0:dt1,:]

# intitialize data frames to save time series data
ser_df_top = pd.DataFrame(index=castdates, columns=data_to_plot)
ser_df_bot = pd.DataFrame(index=castdates, columns=data_to_plot)

# plot the CTD cast data for this station
nr = 0
for cdate in castdates:
    imo = cdate.month
    cast = casts[casts.index==cdate]
    # drop repeat values (aiming for the first of a depth pair)
    zdf = np.diff(cast['Z'])
    zdf = np.concatenate((np.array([1.,]),zdf))
    mask = zdf != 0
    cast = cast[mask]
    cast = cast[:-5] # drop bottom values (sometimes bad)

    for fld in data_to_plot:
        nc = get_nc[fld]
        ax = axes[nr,nc]
        cast.plot(x=fld, y = 'Z', style='-', grid=True,
                  color=month_color_dict[imo], ax=ax, legend=False,
                  xlim=data_range_dict[fld], ylim=data_range_dict['Z'])
        ax.set_xlabel(fld + ' (' + data_unit_dict[fld] + ')')
        
        # also gather a time series entry
        zdiv = -10
        #ser_df_top.loc[cdate, fld] = cast[(cast['Z']>=zdiv) & (cast['Z']<-3)].mean(axis=0)[fld]
        ser_df_top.loc[cdate, fld] = cast[cast['Z']>=zdiv].mean(axis=0)[fld]
        ser_df_bot.loc[cdate, fld] = cast[cast['Z']<zdiv].mean(axis=0)[fld]

# set more things about the cast plots

# add month labels with colors
if nr == 0:
    ax = axes[0, 0]
    for imo in months:
        ax.text(.05, 1 - imo/13.,
            month_name_dict[imo], color=month_color_dict[imo],
            fontsize=14, fontweight='bold',
            verticalalignment='center', transform=ax.transAxes)

fig.suptitle(title)

# NEXT: A figure for TIME SERIES

# make the axes (overwrites previous axes object)
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=figsize,
                         squeeze=False, sharex=True)

# plot data
nr = 0
dft = ser_df_top
dfb = ser_df_bot
for fld in data_to_plot:
    nc = get_nc[fld] 
    ax = axes[nr,nc]
    dft.plot(ax=ax, y=fld, grid=True, style='*-r', legend=False)
    dfb.plot(ax=ax, y=fld, grid=True, style='*-b', legend=False,
             ylim=data_range_dict[fld])
    ax.text(.05, .95, fld + ' (' + data_unit_dict[fld] + ')',
            fontsize=12, fontweight='bold', transform=ax.transAxes)
    ax.set_xlabel('')

# add explanatory text
ax = axes[0, 0]
ax.text(.05, .1, 'Mean above ' + str(zdiv) + ' m', color='r',
    fontsize=14, transform=ax.transAxes)
ax.text(.05, .05, 'Mean below ' + str(zdiv) + ' m', color='b',
    fontsize=14, transform=ax.transAxes)

fig.suptitle(title)

plt.show()
