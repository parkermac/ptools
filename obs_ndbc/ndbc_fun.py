"""
Functions to use for coastal class code.

"""

import pandas as pd
import numpy as np
import os
import shutil

def get_name_dict():
    name_dict = {
         '41002': 'South Hatteras',
         '46002': 'Oregon Offshore',
         '46005': 'Washington Offshore',
         '46014': 'California Shelf',
         '46015': 'Oregon Shelf',
         '46029': 'Washington Shelf',
         '46041': 'Washington Shelf North',
         '46059': 'California Offshore',
         '42040': 'Gulf Coast Mississippi',
         '42035': 'Gulf Coast Texas',
         'sisw1': 'Smith Island',
         'wpow1': 'West Point'}
    return name_dict

def make_dir(dirname, clean=False):
    # Make a directory if it does not exist.
    # Use clean=True to clobber the existing directory.
    if clean == True:
        shutil.rmtree(dirname, ignore_errors=True)
        os.mkdir(dirname)
    else:
        try:
            os.mkdir(dirname)
        except OSError:
            pass # assume OSError was raised because directory already exists

def get_data(fn, tf, yr):
    # read the txt file into a Dataframe
    if yr in range(1970, 1999): # 2-digit year: use %y
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
    if fn == '46002/46002h2015.txt':
        print('** fixing data by hand **')
        df[5800:6250] = np.nan
        
    # create wind stress time series
    # WSPD is in m/s and WDIR is the compass direction
    # that the wind is coming FROM
    
    # First: create 10m standard WSPD
    P = 0.11
    z_stnd = 10
    z_meas = 5
    df['WSPD_10'] = df['WSPD'] * (z_stnd/z_meas)**P
    wspd = df.WSPD_10.values
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