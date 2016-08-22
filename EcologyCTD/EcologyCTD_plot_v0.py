"""
Plots data from a Department of Ecology csv file, for flight CTD stations.

 Data downloaded as .csv from here:
https://fortress.wa.gov/ecy/eap/marinewq/mwdataset.asp
 - at the top choose the station
 - at the bottom choose "csv" and "all years" and then "get file"

10/22/2015 Website moved to:

"""

## USER INPUT ##

# where are the data csv files
# NOTE: Windows users use double forward slashes
# e.g dir0 = 'C://Users//Julie//Downloads//'
dir0 = '/Users/PM5/Documents/tools_data/obs_data/ctd_bottles/DOE_csv/'

# choose which stations to plot by commenting out parts of this list
# names ending in _0 are multi-year 1990-2014
sta_to_plot = [
    #'BLL009_0',
    #'PSB003_0', # edited by hand to remove bad cast 9/18/2002
                # salinity spike 4 m Feb/Mar 1994
                # 2007 many T, s spikes, especially deep
    'DNA001_0', # deep T bad in July 2001
                # 2002 shallow S low value Feb, shallow T spike Oct
                # 2006 lots of little T & s spikes all depths, many months
    'HCB003_0',
    'BUD005_0'
    ]

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

# choose z limit for cast plots
z_deep = -40 # positive up, zero at surface (m)

# set depth limits for averages (time series plot)
ztop = -5 # red curve will be average above this z (m)
zbot = -5 # blue curve will be average below this z (m)

# decide whether or not to print a figure (should appear as a .png
# in whatever directory you are keeping your data files)
make_png = False

# decide if you want to save the casts data to a csv file
do_csv = False

# if this is multi_year data, you can choose to plot a single year
# by making this true, and selecting the year
focus_year = True
fyear = 2013

## END USER INPUT ##

# if desired, print a list of which data files available
if False: # make True to see list of available files
    import os
    flist = os.listdir(dir0)
    print('\nAVAILABLE CSV FILES:\n')
    for fn in flist:
        if 'fixed' not in fn:
            print(fn)

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

# PLOTTING

# setup
import matplotlib.pyplot as plt
plt.close('all')

# plot CASTS

# set up the plot axes (an array)
NR = len(sta_to_plot)
NC = len(data_to_plot)
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(15,8), squeeze=False)

# dictionaries linking stations to row number and data field to column number
get_nr = dict(zip(sta_to_plot, range(NR)))
get_nc = dict(zip(data_to_plot, range(NC)))

# initialize dictionaries for times series data (used in another plot)
ser_top_dict = dict()
ser_bot_dict = dict()

for sta in sta_to_plot:

    # check to see if this is multi-year data
    if sta[-2:] == '_0':
        multi_year = True
    else:
        multi_year = False

    nr = get_nr[sta] # which row to plot in

    # get this station's data frame using our function
    casts = get_casts(sta, dir0, data_name_dict)

    if do_csv:
        casts.to_csv(dir0 + sta + '_processed.csv')

    # identify a single cast by its date
    alldates = casts.index
    castdates = alldates.unique() # a short list of unique dates (1 per cast)
    if multi_year and focus_year:
        from datetime import datetime
        dt0 = datetime(fyear,1,1)
        dt1 = datetime(fyear,12,31)
        castdates = castdates[castdates >= dt0]
        castdates = castdates[castdates <= dt1]

    # intitialize data frames to save time series data
    import pandas as pd
    ser_df_top = pd.DataFrame(index=castdates, columns=data_to_plot)
    ser_df_bot = pd.DataFrame(index=castdates, columns=data_to_plot)

    # plot the CTD cast data for this station
    for cdate in castdates:
        imo = cdate.month
        cast = casts[casts.index==cdate]
        import numpy as np
        # drop repeat values (aiming for the first of a depth pair)
        zdf = np.diff(cast['Z'])
        zdf = np.concatenate((np.array([1.,]),zdf))
        mask = zdf != 0
        cast = cast[mask]
        cast = cast[:-5] # drop bottom values (sometimes bad)

        # here is where I could intervene with a filter like
        #aa = cast['Salinity']
        #aa[(aa - aa.mean()).abs() > 2*aa.std()] = np.nan
        #cast['Salinity'] = aa

        for fld in data_to_plot:
            nc = get_nc[fld]
            axes[nr, nc].plot(cast[fld].values, cast['Z'],'-',
                color=month_color_dict[imo], linewidth = 2)
            # gather a time series entry
            ser_df_top.ix[cdate, fld] = cast[cast['Z']>ztop].mean(axis=0)[fld]
            ser_df_bot.ix[cdate, fld] = cast[cast['Z']<zbot].mean(axis=0)[fld]
    ser_top_dict[sta] = ser_df_top
    ser_bot_dict[sta] = ser_df_bot

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

# you could also save the file to disk with a line like
if make_png:
    plt.savefig(dir0 + 'test_cast.png')

# NEXT: A figure for TIME SERIES

# make the axes (overwrites previous axes object)
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(15,8), squeeze=False)

# plot data
for sta in sta_to_plot:
    nr = get_nr[sta]
    dft = ser_top_dict[sta]
    dfb = ser_bot_dict[sta]
    for fld in data_to_plot:
        nc = get_nc[fld]
        if multi_year and not focus_year:
            axes[nr, nc].plot(dft.index,dft[fld].values,'*-r', linewidth=3)
            axes[nr, nc].plot(dft.index,dfb[fld].values,'o-b', linewidth=3)
        else:
            tt = dft.index.dayofyear
            axes[nr, nc].plot(tt,dft[fld].values,'*-r', linewidth=3)
            axes[nr, nc].plot(tt,dfb[fld].values,'o-b', linewidth=3)

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
        if (not multi_year) or focus_year:
            axes[nr, nc].axis([0, 366, r0, r1])
        else:
            axes[nr, nc].set_ylim(r0, r1)
            from datetime import datetime
            axes[nr, nc].set_xlim(datetime(1990,1,1),datetime(2015,1,1))
        axes[nr, nc].grid() # add grid lines

# add titles and x labels
for fld in data_to_plot:
    nc = get_nc[fld]
    axes[0, nc].set_title(fld + ' (' + data_unit_dict[fld] + ')')
    if multi_year and not focus_year:
        axes[NR - 1, nc].set_xlabel('Date')
    else:
        axes[NR - 1, nc].set_xlabel('Day of year')

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

if make_png:
    plt.savefig(dir0 + 'test_ser.png')

# cleaning up
import os
flist = os.listdir(dir0)
for fn in flist:
    if 'fixed' in fn:
        os.remove(dir0 + fn)