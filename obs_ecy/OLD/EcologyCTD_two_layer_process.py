"""
Process data from a Department of Ecology csv file, for CTD stations.

The main purpose of the code is to create time series of properties
averaged above and below some depth (like 5 m) and save the DataFrames
as pickle files.

Along the way it is a little over-enthusiastic about plotting, but it is
good to see all the casts that go into the averages.

Runs in about a minute for 11 stations.

Data downloaded as .csv from here:
https://fortress.wa.gov/ecy/eap/marinewq/mwdataset.asp
 - at the top choose the station
 - at the BOTTOM choose "csv" and "all years" and then "get file"

"""

import os
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt

## USER INPUT ##

# where are the data csv files
dir0 = '/Users/PM5/Documents/'
indir = dir0 + 'tools_data/obs_data/ctd_bottles/DOE_csv/'
outdir0 = dir0 + 'ptools_output/env_data/'
outdir = outdir0 + 'EcologyCTD/'
try:
    os.mkdir(outdir0)
except OSError:
    pass # assume OSError was raised because directory already exists
try:
    os.mkdir(outdir)
except OSError:
    pass # assume OSError was raised because directory already exists

# choose which stations to plot by commenting out parts of this list
# names ending in _0 are multi-year 1990-2014
sta_to_plot = [
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
# other choices:

# Debugging of input format
if False:
    fn = (indir + sta_to_plot[0] + '.csv')
    ff = open(fn, 'r')
    counter = 0
    for line in ff:
        if counter < 10:
            print(line)
        else:
            break
        counter += 1
    ff.close()

# choose which data fields to plot by commenting out parts of this list
data_to_plot = [
    'Salinity',
    'Temperature',
    'Sigma',
    'Chl',
    'DO',
    #'Trans',
    #'pH'
    ]

# choose z limit for cast plots
z_deep = -25 # positive up, zero at surface (m)

# set depth limits for averages (time series plot)
ztop = -5 # red curve will be average above this z (m)
zbot = -5 # blue curve will be average below this z (m)

## END USER INPUT ##

# lists of data properties

# data long names (used in the csv)
data_long_names = ['salinity(psu)', 'temperature(centigrade)',
    'density(sigmat)', 'chlorophyllraw(ug/l)',
    'dissolvedoxygen(mg/l)adjusted', 'lighttransmission(%)', 'pH']

# data short names, units, and plot ranges
data_names =  ['Salinity','Temperature','Sigma', 'Chl', 'DO',   'Trans',  'pH']
data_units =  ['psu',     'deg C',      'kg/m3', 'ug/l', 'mg/l', '%',     '']
data_ranges = [(14,34),   (4,20),       (14,26), (0,30), (0,16), (0,100), (7,9)]

# dictionaries to look up data attributes using names
data_name_dict = dict(zip(data_long_names, data_names))
data_unit_dict = dict(zip(data_names, data_units))
data_range_dict = dict(zip(data_names, data_ranges))

# more useful lists and dictionaries for plotting different months

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

# a function to read a csv file into a "data frame"
def get_casts(sta, dir0, data_name_dict):
    filename = sta + '.csv'
    filename2 = sta + '_fixed.csv'
    fn = dir0 + filename
    fn2 = dir0 + filename2
    # replace bad data values
    # (hard to parse a csv file is some fields have commas!)
    ff = open(fn, 'r')
    ff2 = open(fn2, 'w')
    for line in ff:
        if '9,999.0' in line:
            line = line.replace('9,999.0','9999.0')
        ff2.write(line)
    ff.close()
    ff2.close()
    # read csv data into a data frame
    import pandas as pd
    casts = pd.read_csv(fn2, parse_dates = [' date'])
    # remove spaces from column headings
    cols = casts.columns
    cols2 = []
    for col in cols:
        Col = col.replace(' ','')
        # and replace column headings with the shorter data names
        if Col in data_name_dict.keys():
            cols2.append(data_name_dict[Col])
        else:
            cols2.append(Col)
    casts.columns = cols2
    # and specify the index column
    casts = casts.set_index('date')
    # and make a Z column
    casts['Z'] = -casts['depth(meters)']
    # and replace missing data with NaN
    import numpy as np
    casts[casts==-9999] = np.nan
    return casts

#%% PLOTTING

fig_size = (15,8)

# setup
plt.close('all')

# plot CASTS

# set up the plot axes (an array)
NR = len(sta_to_plot)
NC = len(data_to_plot)
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(fig_size), squeeze=False)

# dictionaries linking stations to row number and data field to column number
get_nr = dict(zip(sta_to_plot, range(NR)))
get_nc = dict(zip(data_to_plot, range(NC)))

# initialize dictionaries for times series data (used in another plot)
df_top_dict = dict()
df_bot_dict = dict()

for sta in sta_to_plot:

    nr = get_nr[sta] # which row to plot in

    # get this station's data frame using our function
    casts = get_casts(sta, indir, data_name_dict)

    # identify a single cast by its date
    alldates = casts.index
    castdates = alldates.unique() # a short list of unique dates (1 per cast)

    # intitialize data frames to save time series data
    df_top = pd.DataFrame(index=castdates, columns=data_to_plot)
    df_bot = pd.DataFrame(index=castdates, columns=data_to_plot)

    # plot the CTD cast data for this station
    for cdate in castdates:
        imo = cdate.month
        cast = casts[casts.index==cdate]
        # drop repeat values (aiming for the first of a depth pair)
        zdf = np.diff(cast['Z'])
        zdf = np.concatenate((np.array([1.,]),zdf))
        mask = zdf != 0
        cast = cast[mask]
        cast = cast[:-5] # drop bottom values (sometimes bad)[:-5]

        # here is where I could intervene with a filter like
        #aa = cast['Salinity']
        #aa[(aa - aa.mean()).abs() > 2*aa.std()] = np.nan
        #cast['Salinity'] = aa

        for fld in data_to_plot:
            nc = get_nc[fld]
            try:
                axes[nr, nc].plot(cast[fld].values, cast['Z'],'-',
                    color=month_color_dict[imo], linewidth = 2)
                # gather a time series entry
                df_top.ix[cdate, fld] = cast[cast['Z']>=ztop].mean(axis=0)[fld]
                df_bot.ix[cdate, fld] = cast[cast['Z']<zbot].mean(axis=0)[fld]
            except:
                continue

    df_top_dict[sta] = df_top
    df_bot_dict[sta] = df_bot

# set more things about the cast plots

for sta in sta_to_plot:
    nr = get_nr[sta]

    # add axes labels
    axes[nr, 0].set_ylabel('Z (m)')
    if nr == NR -1:
        for fld in data_to_plot:
            nc = get_nc[fld]
            axes[nr, nc].set_xlabel(fld + ' (' + data_unit_dict[fld] + ')')

    # set axes limits, same for each data type
    for fld in data_to_plot:
        nc = get_nc[fld]
        r0, r1 = data_range_dict[fld]
        axes[nr, nc].axis([r0, r1, z_deep, 0])

    # add month labels with colors
    if nr == 0:
        ax = axes[0, 0]
        for imo in months:
            ax.text(.05, 1 - imo/13.,
                month_name_dict[imo], color=month_color_dict[imo], fontsize=10,
                verticalalignment='center', transform=ax.transAxes)

    # add station name
    ax = axes[nr, NC-1]
    ax.text(.95, .05, sta, fontsize=14,
        horizontalalignment='right',
        transform=ax.transAxes)

# then make the plot appear on the screen
plt.show()

#%% NEXT: A figure for TIME SERIES of top and bottom LAYERS

# make the axes (overwrites previous axes object)
NR = len(sta_to_plot)
NC = len(data_to_plot)
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(fig_size), squeeze=False)

# plot data
for sta in sta_to_plot:
    nr = get_nr[sta]
    dft = df_top_dict[sta]
    dfb = df_bot_dict[sta]
    for fld in data_to_plot:
        nc = get_nc[fld]
        axes[nr, nc].plot(dft.index,dft[fld].values,'o-r', linewidth=1)
        axes[nr, nc].plot(dft.index,dfb[fld].values,'o-b', linewidth=1)

# add station names
for sta in sta_to_plot:
    nr = get_nr[sta]
    ax = axes[nr, NC-1]
    ax.text(.95, .95, sta, fontsize=14,
        horizontalalignment='right',
        verticalalignment='top',
        transform=ax.transAxes)
    # set axes limits
    for fld in data_to_plot:
        nc = get_nc[fld]
        r0, r1 = data_range_dict[fld]
        axes[nr, nc].set_ylim(r0, r1)
        axes[nr, nc].set_xlim(datetime(1990,1,1),datetime(2015,1,1))
        axes[nr, nc].grid() # add grid lines

# add titles and x labels
for fld in data_to_plot:
    nc = get_nc[fld]
    axes[0, nc].set_title(fld + ' (' + data_unit_dict[fld] + ')')
    axes[NR - 1, nc].set_xlabel('Date')

# add explanatory text
ax = axes[0, 0]
ax.text(.05, .2,
    'Mean above ' + str(ztop) + ' m', color='r',
    fontsize=14,
    transform=ax.transAxes)
ax.text(.05, .05,
    'Mean below ' + str(zbot) + ' m', color='b',
    fontsize=14,
    transform=ax.transAxes)

plt.show()

#%% NEXT: A figure for TIME SERIES of DENSITY DIFFERENCE by SEASON
NC = 1
# make the axes (overwrites previous axes object)
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(fig_size), squeeze=False)

# plot data
for sta in sta_to_plot:
    nr = get_nr[sta]
    dft = df_top_dict[sta]
    dfb = df_bot_dict[sta]
    delta = dfb - dft
    years = range(1990, 2015)
    #for fld in ['Sigma']:
    #    axes[nr, NC-1].plot(delta.index,delta[fld].values,'-k', linewidth=2, alpha=0.2)
    # add half-year averages
    year_ax = []
    dmean = []
    for year in years:
        dd0 = delta['Sigma'][delta.index.year == year - 1]
        dd1 = delta['Sigma'][delta.index.year == year]
        dd00 = dd0[dd0.index.month > 9]
        dd11 = dd1[dd1.index.month < 4]
        ddd = pd.concat([dd00, dd11])
        dmean.append(ddd.mean())
        year_ax.append(datetime(year,1,1))
    axes[nr, NC-1].plot(year_ax,dmean,'-ob', linewidth=3)
    year_ax = []
    dmean = []
    for year in years:
        dd0 = delta['Sigma'][delta.index.year == year]
        dd00 = dd0[dd0.index.month >=4]
        ddd = dd00[dd00.index.month <= 9]
        dmean.append(ddd.mean())
        year_ax.append(datetime(year,6,30))
    axes[nr, NC-1].plot(year_ax,dmean,'-or', linewidth=3)

# add station names
for sta in sta_to_plot:
    nr = get_nr[sta]
    ax = axes[nr, NC-1]
    ax.text(.95, .95, sta, fontsize=14,
        horizontalalignment='right',
        verticalalignment='top',
        transform=ax.transAxes)
    # set axes limits
    for fld in data_to_plot:
        nc = NC-1
        #axes[nr, nc].set_ylim(0, 8)
        axes[nr, nc].set_xlim(datetime(1990,1,1),datetime(2015,1,1))
        axes[nr, nc].grid() # add grid lines

fld = 'Sigma'
nc = NC-1
axes[0, nc].set_title(fld + ' (' + data_unit_dict[fld] + ')')
axes[NR - 1, nc].set_xlabel('Date')

# add explanatory text
ax = axes[0, 0]
ax.text(.05, .95,
    '(Mean below ' + str(zbot) + ' m) - (Mean above ' + str(ztop) + ' m)', color='k',
    fontsize=14,
    transform=ax.transAxes)
ax.text(.07, .75, 'Summer', color='r', fontsize=16, transform=ax.transAxes)
ax.text(.07, .7, 'Winter', color='b', fontsize=16, transform=ax.transAxes)

plt.show()

# save some results
for sta in sta_to_plot:
    dft = df_top_dict[sta]
    dfb = df_bot_dict[sta]
    dft.to_pickle(outdir + 'df_top_' + sta + '.p')
    dfb.to_pickle(outdir + 'df_bot_' + sta + '.p')

# cleaning up
flist = os.listdir(indir)
for fn in flist:
    if 'fixed' in fn:
        os.remove(indir + fn)