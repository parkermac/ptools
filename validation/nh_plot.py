# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 16:18:48 2016

@author: PM5

Code to read in NH mooring data and compare with a mooring extraction.

"""

#%% setup
import os
import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()
import zfun

import pandas as pd
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import os
import netCDF4 as nc4
import numpy as np

#%% gather mooring data

whichmoor = 'NH10' # 'NH10' or 'NH20'

indir = (Ldir['parent'] + '/tools_data/obs_data/mooring/NH_Line/'
            + whichmoor +'/')
a = os.listdir(indir)

# ********** CTD **********
# header lines:
# DOY_SBE16 Press_SBE16 Temp_SBE16 Sal_SBE16 O2_opt iyr
# DOY_SBE16 Press_SBE16 Temp_SBE16 Sal_SBE16 O2_opt SBE16Yr
fn_list = []
for item in a:
    if 'SBE16' in item:
        fn_list.append(item)
try:
    fn_list.remove('NH20_BBL_SBE16_2013') # has no header
except ValueError:
    pass
DF = pd.DataFrame()
for fn in fn_list:
    print('\n' + fn)
    count = 0
    for line in open(indir + fn, errors='ignore'):
        if 'DOY_SBE16' in line:
            n_header = count
        if count < 20:
            print(line, end='')
        count += 1
    df = pd.read_csv(indir + fn, delim_whitespace=True, header=n_header)
    if 'SBE16Yr' in df.columns:
        yr_name = 'SBE16Yr'
    elif 'iyr' in df.columns:
        yr_name = 'iyr'
    dt = []
    for ii in df.index:
        dt.append( datetime(df.ix[ii,yr_name],1,1) +
            timedelta(df.ix[ii,'DOY_SBE16']-1) )
    df.index = dt
    DF = pd.concat([DF, df])
DF = DF.sort_index()
DFsbe = DF

# ********** C02 **********
# headers lines:
# Day_of_Year     SAMI_Temp       SAMI_CO2        SAMI_DateTime   SAMIYr
fn_list = []
for item in a:
    if 'CO2' in item:
        fn_list.append(item)
DF = pd.DataFrame()
for fn in fn_list:
    print('\n' + fn)
    count = 0
    for line in open(indir + fn, errors='ignore'):
        if 'Day_of_Year' in line:
            n_header = count
        if count < 20:
            print(line, end='')
        count += 1
    cols=['Day_of_Year','SAMI_Temp','SAMI_CO2','SAMI_Date','Time','SAMIYr']
    df = pd.read_csv(indir + fn, delim_whitespace=True,
                     skiprows=n_header+1, names=cols)
    if 'SAMIYr' in df.columns:
        yr_name = 'SAMIYr'
    dt = []
    for ii in df.index:
        dt.append( datetime(df.ix[ii,yr_name],1,1) +
            timedelta(df.ix[ii,'Day_of_Year']-1) )
    df.index = dt
    DF = pd.concat([DF, df])
try:
    DF.ix[:,'SAMI_CO2'][DF.ix[:,'SAMI_CO2'] == 'NaN4'] = np.nan
except TypeError:
    pass
DF = DF.dropna()
DF.ix[:,'SAMI_CO2'] = DF.ix[:,'SAMI_CO2'].astype(float)
DF = DF.sort_index()
DFco2 = DF

# ********** pH **********
# header lines:
# Day_of_Year     SAMI_Temp       SAMI_pH SAMI_DateTime   SAMI_year
fn_list = []
for item in a:
    if 'pH' in item:
        fn_list.append(item)
DF = pd.DataFrame()
for fn in fn_list:
    print('\n' + fn)
    count = 0
    for line in open(indir + fn, errors='ignore'):
        if 'Day_of_Year' in line:
            n_header = count
        if count < 20:
            print(line, end='')
        count += 1
    cols=['Day_of_Year','SAMI_Temp','SAMI_pH','SAMI_Date','Time','SAMIYr']
    df = pd.read_csv(indir + fn, delim_whitespace=True,
                     skiprows=n_header+1, names=cols)
    if 'SAMIYr' in df.columns:
        yr_name = 'SAMIYr'
    elif 'SAMI_yr' in df.columns:
        yr_name = 'SAMI_yr'
    dt = []
    for ii in df.index:
        dt.append( datetime(df.ix[ii,yr_name],1,1) +
            timedelta(df.ix[ii,'Day_of_Year']-1) )
    df.index = dt
    DF = pd.concat([DF, df])
DF = DF.sort_index()
DFph = DF

#%% get the mooring extraction

moor_dir = Ldir['LOo'] + 'moor/'
moor_fn = ('cascadia1_base_lobio1_'+ whichmoor +
            '_low_pass_2013.01.02_2015.12.30.nc')
ds = nc4.Dataset(moor_dir + moor_fn)
z = ds['z_rho'][:,0]
zz = np.min(z) + 2
izz = zfun.find_nearest_ind(z, zz)
temp = ds['temp'][izz, :]
salt = ds['salt'][izz, :]
oxygen = ds['oxygen'][izz, :]
pH = ds['PH'][izz, :]
ot = ds['ocean_time'][:]
ds.close()
mdt = []
for item in ot:
    mdt.append(datetime(1970,1,1) + timedelta(days=item/86400))
DFmoor = pd.DataFrame(index=mdt, columns=['Temp', 'Salt', 'Oxygen', 'pH'])
DFmoor['temp'] = temp
DFmoor['salt'] = salt
DFmoor['oxygen'] = oxygen
DFmoor['pH'] = pH

#%% plotting
plt.close('all')
fig = plt.figure(figsize=(18,10))
mks = 2

ax = fig.add_subplot(221)
DFsbe.plot(y='Temp_SBE16', style='.k', ax=ax, grid=True, markersize=mks)
DFco2.plot(y='SAMI_Temp', style='.r', ax=ax, grid=True, markersize=mks)
DFph.plot(y='SAMI_Temp', style='.g', ax=ax, grid=True, markersize=mks)
DFmoor.plot(y='temp', style='-b', ax=ax, grid=True)
ax.legend(['sbe','co2','ph'])
ax.set_ylabel('Temperature (C)')
ax.set_title(indir)
xl = ax.get_xlim()

ax = fig.add_subplot(222)
DFsbe.plot(y='O2_opt', style='.k', ax=ax, grid=True, markersize=mks)
DFmoor.plot(y='oxygen', style='-b', ax=ax, grid=True)
ax.legend(['sbe','roms'])
ax.set_ylabel('Oxygen')
ax.set_xlim(xl)

ax = fig.add_subplot(223)
DFsbe.plot(y='Sal_SBE16', style='.k', ax=ax, grid=True, markersize=mks)
DFmoor.plot(y='salt', style='-b', ax=ax, grid=True)
ax.legend(['sbe','roms'])
ax.set_ylabel('Salinity')
ax.set_xlim(xl)

ax = fig.add_subplot(224)
DFph.plot(y='SAMI_pH', style='.g', ax=ax, grid=True, markersize=mks)
DFmoor.plot(y='pH', style='-b', ax=ax, grid=True)
ax.legend(['sami','roms'])
ax.set_ylabel('pH')
ax.set_xlim(xl)



