"""
Code to plot observed tide time series.

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np

from importlib import reload
import transit
reload(transit)

def read_tide(in_fn, city):
    df = pd.read_csv(in_fn, index_col='Date Time', parse_dates = True)
    for k in df.keys():
        df = df.rename(columns={k: k.strip()})
    df = df.drop(['Sigma', 'I', 'L'], axis=1)
    df = df.rename(columns={'Water Level': city})
    return df


home = os.environ.get('HOME')
dir00 = home + '/Documents/'

# read in tide data (time is UTC)
indir = dir00 + 'ptools_data/tide/'

fn = 'CO-OPS__9447130__hr.csv'
city = 'Seattle'
df1 = read_tide(indir+fn, city)

fn = 'CO-OPS__9441102__hr.csv'
city = 'Westport'
df2 = read_tide(indir+fn, city)
    

# PLOTTING LOOP

# initializing
plt.close('all')

fig = plt.figure(figsize=(14,8))
ax = fig.add_subplot(111)
lw = .5
df1.plot(ax=ax, lw=lw)
df2.plot(ax=ax, lw=lw)

plt.show()
