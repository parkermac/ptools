# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 13:10:34 2016

@author: PM5

Code to plot river tracks and names, and then associate the names with
rivers that have temperature climatology data.

"""

import os
import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)

from importlib import reload
import Lfun
reload(Lfun)
Ldir = Lfun.Lstart(gridname='test')

import gfun; reload(gfun)
G = gfun.gstart()
import pfun

import pandas as pd
import matplotlib.pyplot as plt

rivp = os.path.abspath('../../LiveOcean/forcing/riv1/')
if rivp not in sys.path:
    sys.path.append(rivp)
import river_functions as rivfun
reload(rivfun)

#%% get river info
ri_fn = G['ri_dir'] + 'river_info.csv'

df = pd.read_csv(ri_fn, index_col='rname')

dir0 = Ldir['data'] + 'rivers/'
clim_dir = dir0 + 'Data_T_clim/'

c_list_raw = os.listdir(clim_dir)
c_list = []
for m in c_list_raw:
    if '.csv' in m:
        c_list.append(m)

#%% associate rivers with ones that have temperature climatology data

df = rivfun.get_tc_rn(df)

#%% plotting
plt.close()

fig = plt.figure(figsize=(12,17))
ax = fig.add_subplot(111)
ax_lims = [-127.5, -121, 42.5, 50.5]
ax.set_xlim(ax_lims[:2])
ax.set_ylim(ax_lims[-2:])
pfun.add_coast(ax)
pfun.dar(ax)

for rn in df.index:
    try:
        fn_tr = G['ri_dir'] + 'tracks/' + rn + '.csv'
        df_tr = pd.read_csv(fn_tr, index_col='ind')
        x = df_tr['lon'].values
        y = df_tr['lat'].values
        ax.plot(x, y, '-r', linewidth=2)
        ax.plot(x[-1], y[-1], '*r')

        if (rn + '.csv') in c_list:
            ax.text(x[-1]+.05, y[-1]+.03, rn.title(), color='r', weight='bold')
        else:
            ax.text(x[-1]+.05, y[-1]+.03, rn.title())
    except FileNotFoundError:
        pass

ax.set_title('Rivers with Flow Data (RED have Temp Climatology)')

plt.show()