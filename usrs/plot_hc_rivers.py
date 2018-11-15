#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot as-run river time series.

"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime

import os
import sys

pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import zrfun
import zfun

pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

pth = os.path.abspath('../pgrid')
if pth not in sys.path:
    sys.path.append(pth)
import gfun

pth = os.path.abspath('../../LiveOcean/forcing/riv1/')
if pth not in sys.path:
    sys.path.append(pth)
import river_functions as rivfun

# get river flow for a specific year
Ldir = Lfun.Lstart('cas4', 'v2')
fnr = 'cas4_v2_2017.01.01_2017.12.31.p'
fn = Ldir['LOo'] + 'river/' + fnr
df = pd.read_pickle(fn)
# set limits for plotting
rind = df.index
dt0 = rind[0]
dt1 = rind[-1]

# times to emphasize
dt00 = rind[90]
dt11 = rind[90+30]

# get river flow climatology
df_clim = pd.DataFrame(index=df.index, columns=df.columns)
for riv_name in df.columns:
    clm_fn = Ldir['data'] + 'rivers/Data_clim/' + riv_name + '.csv'
    dfc = pd.read_csv(clm_fn, header=None, index_col=0, names=['Qr'])
    df_clim[riv_name] = dfc.loc[:,'Qr'].values

# get river paths for map
G = gfun.gstart()


#%% get river info
ri_fn = G['ri_dir'] + 'river_info.csv'
df_map = pd.read_csv(ri_fn, index_col='rname')

# plotting
plt.close('all')
fig = plt.figure(figsize=(14,8))

rn_list = ['dosewallips', 'duckabush', 'hamma', 'skokomish']
rn_dict = dict(zip(rn_list,[1,2,4,5]))

# river flow time series plots
rr = 1
q0 =0
q1 = 150
for riv_name in rn_list:
    ax = fig.add_subplot(2,3,rn_dict[riv_name])
    if rr==3:
        do_legend = True
    else:
        do_legend = False
    df.loc[:,riv_name].plot(ax=ax, color='cornflowerblue',
        label='Observed 2017', legend=do_legend)
    df_clim.loc[:,riv_name].plot(ax=ax, color='orange',
        label='Climatology', legend=do_legend)
    df.loc[dt00:dt11,riv_name].plot(ax=ax, color='b', linewidth=3,
        label='Obs. April', legend=do_legend)
    df_clim.loc[dt00:dt11,riv_name].plot(ax=ax, color='r',
        label='Clim. April', linewidth=3, legend=do_legend)
    ax.set_xlim(dt0, dt1)
    ax.set_ylim(q0, q1)
    ax.text(.05, .9, riv_name.title(), transform=ax.transAxes,
        horizontalalignment='left', fontweight='bold')
    #ax.set_title(riv_name.title())
    if rr in [1,3]:
        ax.set_ylabel('Flow $(m^{3}\ s^{-1})$')
    if rr in [1,2]:
        ax.set_xticklabels([])
    else:
        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
        ax.xaxis.set_tick_params(labelrotation=45)
        ax.set_xlabel('Date 2017')
    # add emphasis for selected month
    rr += 1
    
# river location map
ax = fig.add_subplot(1,3,3)
aa = [-123.2, -122.6, 47, 48]
ax.axis(aa)
pfun.add_coast(ax)
pfun.dar(ax)

for rn in rn_list:
    fn_tr = G['ri_dir'] + 'tracks/' + rn + '.csv'
    df_tr = pd.read_csv(fn_tr, index_col='ind')
    x = df_tr['lon'].values
    y = df_tr['lat'].values
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[-1], y[-1], '*r')
    ax.text(x[-1], y[-1], rn.title(), fontweight='bold')
ax.set_title('Hood Canal Rivers')

plt.show()
