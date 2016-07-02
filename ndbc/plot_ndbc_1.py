# -*- coding: utf-8 -*-
"""
Plot data on river flow and nearby NDBC buoy wave data.

"""

# import modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

# **** USER EDITS ****

sn = 'wpow1'


# set time limits to plot
if True: # plot the full record of data I downloaded
    t0 = datetime(1984,1,1)
    t1 = datetime(2015,12,31)
else: # or select a smaller time span, like a year
    t0 = datetime(2011,1,1)
    t1 = datetime(2012,12,31)

# choose time filtering length (a string) e.g.:
#    m = month
#    w = week
#    d = day
#    h = hour (no filtering)
tf = 'm'
# but for this exercise it is probably fine to leave this at 'd'

# **** END USER EDITS ****


# load and process NDBC Buoy Data

"""
For wtpo1 the columns are:

1984-1993:
YY MM DD hh WD   WSPD GST  WVHT  DPD   APD  MWD  BAR    ATMP  WTMP  DEWP  VIS

1994-1995:
YY MM DD hh WD   WSPD GST  WVHT  DPD   APD  MWD  BAR    ATMP  WTMP  DEWP  VIS\r

1996-1997:
YY MM DD hh WD   WSPD GST  WVHT  DPD   APD  MWD  BAR    ATMP  WTMP  DEWP  VIS

1998:
YY MM DD hh WD   WSPD GST  WVHT  DPD   APD  MWD  BAR    ATMP  WTMP  DEWP  VIS
98 01 01 00 172  8.2  8.9 99.00 99.00 99.00 999 1011.4   8.7 999.0   6.9 99.0\r

1999:
YYYY MM DD hh WD   WSPD GST  WVHT  DPD   APD  MWD  BAR    ATMP  WTMP  DEWP  VIS
1999 01 01 00  69  4.5  6.1 99.00 99.00 99.00 999 1021.6   7.6 999.0   4.8 99.0\r

2000-2004:
YYYY MM DD hh WD   WSPD GST  WVHT  DPD   APD  MWD  BAR    ATMP  WTMP  DEWP  VIS  TIDE

2005-2006:
YYYY MM DD hh mm  WD  WSPD GST  WVHT  DPD   APD  MWD  BAR    ATMP  WTMP  DEWP  VIS  TIDE

2007-2015: (but 2008 and 2010 had deg instead of degT, just for MWD)
#YY  MM DD hh mm WDIR WSPD GST  WVHT   DPD   APD MWD   PRES  ATMP  WTMP  DEWP  VIS  TIDE
#yr  mo dy hr mn degT m/s  m/s     m   sec   sec degT   hPa  degC  degC  degC   mi    ft
"""

def get_data(fn, tf, yr):
    # read the txt file into a Dataframe

    if yr in range(1984, 1999): # 2-digit year: use %y
        df = pd.read_csv(fn, delim_whitespace=True, index_col='date',
                     parse_dates={'date':[0, 1, 2, 3]},
                     date_parser=lambda x: pd.datetime.strptime(x,
                     '%y %m %d %H'))
    if yr in range(1999, 2005):
        df = pd.read_csv(fn, delim_whitespace=True, index_col='date',
                     parse_dates={'date':[0, 1, 2, 3]},
                     date_parser=lambda x: pd.datetime.strptime(x,
                     '%Y %m %d %H'))
    elif yr >= 2005: # add minutes column
        df = pd.read_csv(fn, delim_whitespace=True, index_col='date',
                     skiprows=[1],
                     parse_dates={'date':[0, 1, 2, 3, 4]},
                     date_parser=lambda x: pd.datetime.strptime(x,
                     '%Y %m %d %H %M'))
    df = df.rename(columns={'WD': 'WDIR', 'BAR': 'PRES'})


    # mask known missing data
    df[df==9999.0] = np.nan
    df[df==999.0] = np.nan
    df[df==99.0] = np.nan
    # fix some obviously bad data
    if fn == '46002h2015.txt':
        df[5800:6250] = np.nan
    # create wind stress time series
    # WSPD is in m/s and WDIR is the compass direction
    # that the wind is coming FROM
    # ISSUE: this may not be 10 m wind, so Cd is incorrect
    wspd = df.WSPD.values
    wdir = df.WDIR.values
    theta = 1.5*np.pi - np.pi*wdir/180.
    Cd = 0.0013
    rho_air = 1.22
    tau = Cd * rho_air * wspd**2
    taux = tau * np.cos(theta)
    tauy = tau * np.sin(theta)
    df['taux'] = taux
    df['tauy'] = tauy
    if tf == 'm':
        loff = '-15d'
    elif tf == 'w':
        loff = '-3d'
    else:
        loff = None
    if tf == 'h':
        dff = df
    else:
        dff = df.resample(tf, how='mean', loffset=loff)
    dff = dff[dff.index.year == yr]
    return dff

# process all the NDBC data
yr_list = range(1984,2016)
id_list = []
count = 0
for yr in yr_list:
    id = str(sn) + 'h' + str(yr)
    id_list.append(id)
    fn = ('/Users/PM5/Documents/tools_data/obs_data/ndbc/'
          + sn + '/' + id + '.txt')
    try:
        dff = get_data(fn, tf, yr)
        print(fn + ' = success')
        if count == 0:
            DFF = dff
            count += 1
        else:
            DFF = DFF.append(dff)
    except:
        print(fn + ' = fail')
        pass

# create a dict relating varible names to their units
header = pd.read_csv(fn, nrows=1, delim_whitespace=True)
unit_dict = dict(header.ix[0])
unit_dict['taux'] = 'Pa'
unit_dict['tauy'] = 'Pa'
unit_dict['wspd_clim'] = unit_dict['WSPD']

#%% make climatology

DFF['wspd_clim'] = np.nan
for t in DFF.index:
    mo = t.month
    the_mean = DFF.ix[DFF.index.month==mo, 'WSPD'].mean(axis=0)
    DFF.ix[t, 'wspd_clim'] = the_mean

#%% PLOTTING

plt.close()

# (1) time series figure

fig1 = plt.figure(figsize = (14, 10))
NC = 3
ax1 = fig1.add_subplot(NC,1,1)
ax2 = fig1.add_subplot(NC,1,2)
ax3 = fig1.add_subplot(NC,1,3)


ax = ax1
vn = 'WSPD'
DFF[vn].plot(ax=ax, grid=True, xlim=(t0,t1), ylim=(0, 10), color='r')
vn2 = 'wspd_clim'
DFF[vn2].plot(ax=ax, grid=True, xlim=(t0,t1), ylim=(0, 10), color='b')
ax.set_xticklabels([])
ax.set_xlabel('')
ax.set_ylabel(vn + ' [' + unit_dict[vn] + ']')
ax.text(.05, .85, 'NDBC Station ' + str(sn), transform=ax1.transAxes)

ax = ax2
vn = 'tauy'
DFF[vn].plot(ax=ax, grid=True, xlim=(t0,t1), ylim=(-.05, .15))
ax.set_xticklabels([])
ax.set_xlabel('')
ax.set_ylabel(vn + ' [' + unit_dict[vn] + ']')

ax = ax3
vn = 'ATMP'
DFF[vn].plot(ax=ax, grid=True, xlim=(t0,t1), ylim=(0, 20))
ax.set_ylabel(vn + ' [' + unit_dict[vn] + ']')

plt.show()
