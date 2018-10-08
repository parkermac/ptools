#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots data from the Department of Ecology, combining data from CTD casts
and bottles at the same station.  Also including model output.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()
import zfun

Ldir['gtagex'] = 'cas4_v2_lo6biom'
year = 2017
dir11 = Ldir['parent'] + 'ptools_output/ecology/'
out_fn = 'ObsMod_' + Ldir['gtagex'] + '_'+ str(year) + '.p'
print('Loading ' + out_fn)
Bc = pd.read_pickle(dir11 + out_fn)

lim_dict = {'Salinity': (0, 34), 'Temp. (deg C)': (0, 24),
    'DO (mg L-1)': (0, 20), 'DIN (uM)': (0,55), 'Chl (mg m-3)': (0,50)}
Z_list = [0,-10,-30]
clist = ['r', 'g', 'b']

# Print some statistics to the screen
def errstat(Bc, vn):
    ObsMod = Bc[[vn, 'Mod '+vn]]
    ObsMod = ObsMod.dropna()
    count = len(ObsMod)
    rmse = np.sqrt( ((Bc['Mod '+vn] - Bc[vn])**2).mean() )
    bias = (Bc['Mod '+vn] - Bc[vn]).mean()
    errstr = ('RMSE = %5.2f, Bias = %5.2f (N=%4d)' % (rmse, bias, count))
    return errstr

plt.close('all')
# make overall scatterplots
fig = plt.figure(figsize=(14,9))
NR = 2; NC = 3
pp = 1
for vn in ['Salinity', 'Temp. (deg C)', 'DO (mg L-1)', 'DIN (uM)', 'Chl (mg m-3)']:
    ax = fig.add_subplot(NR, NC ,pp)
    ii = 0
    for Zn in Z_list:
        Bcz = Bc[Bc['Znom']==Zn]
        Bcz.plot(x=vn, y = 'Mod '+vn, style = 'o'+ clist[ii],
            grid=True, ax=ax, legend=False, alpha=.4)
        ii += 1
    ax.set_xlim(lim_dict[vn])
    ax.set_ylim(lim_dict[vn])
    ax.set_title(errstat(Bc, vn))
    v0 = lim_dict[vn][0]; v1 = lim_dict[vn][1]
    ax.plot([v0, v1], [v0, v1], '-k')
    ax.text(.05, .9, vn, fontweight='bold', transform=ax.transAxes)
    ax.set_xlabel('Observations')
    ax.set_ylabel('Model')
    
    if pp==1:
        ax.text(.05, .5, 'Z =   0 m', color='r', fontweight='bold', transform=ax.transAxes)
        ax.text(.05, .4, 'Z = -10 m', color='g', fontweight='bold', transform=ax.transAxes)
        ax.text(.05, .3, 'Z = -30 m', color='b', fontweight='bold', transform=ax.transAxes)
    
    pp+=1
plt.show()