"""
Plot a weather record from NCDC.
"""

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np

# SSMSP import
import os
import sys
pth = os.path.abspath('../ssmsp')
if pth not in sys.path:
    sys.path.append(pth)
import sfun
from importlib import reload
reload(sfun)

in_dir = '../../ptools_data/ncdc/'
fn = '1527887.csv' # SeaTac record

out_dir = '../../ptools_output/ncdc/'
sfun.make_dir(out_dir)

a = pd.read_csv(in_dir + fn, low_memory=False)

a['DATE'] = pd.to_datetime(a['DATE'])
a = a.set_index('DATE')

# form monthly averages
am = a.resample('M').mean()

# Processing rain for phenology
# wateryear time: starts on October 1 of the year before, so wyt=1/1/1949
# happens on DATE=10/1/1948
wyt = a.index + timedelta(days=92)
prcp = pd.DataFrame(a['PRCP'].copy())
prcp['wyt'] = wyt
prcp = prcp.set_index('wyt')

# make the sum over each wateryear, to highlight wet and dry years
net_pr = prcp.resample('Y').sum()
net_pr = net_pr.iloc[1:-2] # drop first and last two years
net_pr = net_pr/10 # convert mm to cm
net_pr.index = net_pr.index.year

# now make a dataframe where the index is water year yearday
# and each column is the cumsum of PRCP for that water year
wy0 = 1949
wy1 = 2017

cs_pr = pd.DataFrame(index=range(1,366))
cs_pr.index.name = 'Yearday'
for wy in range(wy0, wy1+1):
    this_ser = prcp[prcp.index.year == wy]
    this_ser.index = this_ser.index.dayofyear
    cs_pr[wy] = this_ser.cumsum()/1000

# form time series of cumulative rain normalized by the
# total for that year
cs_pr_sc = cs_pr.copy()/cs_pr.loc[365,:]

# then pick out the yearday of certain scaled net rain:
dd = pd.DataFrame(index=range(wy0, wy1+1))
dd['d15'] = (cs_pr_sc <= .15).sum()
dd['d85'] = (cs_pr_sc <= .85).sum()

# PLOTTING
fs=16
plt.close('all')

# PRECIPITATION PLOT
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)

dd['d85'].plot(ax=ax, label='Yearday for 85% of Total',
    color='darkorange', linewidth=2, legend=True)
    
net_pr['PRCP'].plot(ax=ax, label='Net Precipitation (cm)',
    style='-b', linewidth=4, alpha=.5, legend=True)
    
dd['d15'].plot(ax=ax, label='Yearday for 15% of Total',
    color='purple', linewidth=2, legend=True, grid=True)
    
ax.set_xlabel('Water Year (Oct.-Sept.)', fontsize=fs)

ax.tick_params(labelsize=fs-2)
plt.setp(plt.gca().get_legend().get_texts(), fontsize=fs)

ax.set_xlim(wy0, wy1)
ax.set_title('SeaTac Precipitation Data', fontsize=fs+2)

plt.savefig(out_dir + 'Seatac_rain.png')

# TEMPERATURE and SUNSHINE PLOT
fig = plt.figure(figsize=(12,10))

# sunshine
ax = fig.add_subplot(211)
am['TSUN'] = am['TSUN']/60 # convert minutes of sunshine per day to hours
am['TSUN'].plot(ax=ax, label='Hours of Sunshine per day', grid=True,
    legend=True, color='red')
ax.set_xlabel('')
#ax.set_xticklabels([])
ax.set_xlim(datetime(wy0,1,1), datetime(wy1,12,31))
ax.set_title('SeaTac Monthly Sunshine and Temperature Extremes', fontsize=fs+2)
ax.tick_params(labelsize=fs-2)
plt.setp(plt.gca().get_legend().get_texts(), fontsize=fs)

# temperatures
ax = fig.add_subplot(212)
am['TMAX'].plot(ax=ax, grid=True, label='Daily Max Temperature (degC)',
    legend=True, color='darkorange')
am['TMIN'].plot(ax=ax, grid=True, label='Daily Min Temperature (degC)',
    legend=True, color='blue')
ax.tick_params(labelsize=fs-2)
plt.setp(plt.gca().get_legend().get_texts(), fontsize=fs)
ax.set_xlim(datetime(wy0,1,1), datetime(wy1,12,31))
ax.set_xlabel('Date', fontsize=fs)

plt.show()

plt.savefig(out_dir + 'Seatac_Temperature_Sunshine.png')
