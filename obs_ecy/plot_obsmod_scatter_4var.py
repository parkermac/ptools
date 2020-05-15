#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots data from the Department of Ecology, combining data from CTD casts
and bottles at the same station.  Also including model output.

Just like plot_obsmod_scatter, but only plots 4 variables instead of 6.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seawater as sw

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()
import zfun

Ldir['gtagex'] = 'cas6_v3_lo8b'
year = 2017
dir11 = Ldir['parent'] + 'ptools_output/ecology/'
out_fn = 'ObsMod_' + Ldir['gtagex'] + '_'+ str(year) + '.p'
print('Loading ' + out_fn)
Bc = pd.read_pickle(dir11 + out_fn)

# output plot name
out_plot_fn = dir11 + 'ecy_scatter_4var_' + Ldir['gtagex'] + '_'+ str(year) + '.png'

lim_dict = {'Salinity': (15,35), 'Temp. (deg C)': (5,20),
    'DO (mg L-1)': (0,12), 'DIN (uM)': (0,35)}
Z_list = [0,-10,-30, -80]
clist = ['r', 'g', 'b', 'k']

# Print some statistics to the screen
def errstat(Bc, vn):
    ObsMod = Bc[[vn, 'Mod '+vn]]
    ObsMod = ObsMod.dropna()
    count = len(ObsMod)
    rmse = np.sqrt( ((Bc['Mod '+vn] - Bc[vn])**2).mean() )
    bias = (Bc['Mod '+vn] - Bc[vn]).mean()
    errstr = ('RMSE = %5.2f\nBias = %5.2f\n(N=%4d)' % (rmse, bias, count))
    return errstr

# PLOTTING
fs = 18
abc = 'abcdef'
balpha=.7
vn_dict = {'Salinity':'Salinity', 'Temp. (deg C)':'Temperature [$^{\circ} C$]',
    'DO (mg L-1)':'DO [$mg\ L^{-1}$]', 'DIN (uM)':'DIN [$\mu M$]'}
plt.close('all')
# make overall scatterplots
fig = plt.figure(figsize=(12,12))
NR = 2; NC = 2
pp = 1
for vn in ['Salinity', 'Temp. (deg C)', 'DO (mg L-1)', 'DIN (uM)']:
    ax = fig.add_subplot(NR, NC ,pp)
    ii = 0
    for Zn in Z_list:
        Bcz = Bc[Bc['Znom']==Zn]
        try:
            Bcz.plot(x=vn, y = 'Mod '+vn, style = 'o'+ clist[ii],
                grid=True, ax=ax, legend=False, alpha=.4)
        except ValueError:
            pass
        ii += 1
    ax.set_xlim(lim_dict[vn])
    ax.set_ylim(lim_dict[vn])
    v0 = lim_dict[vn][0]; v1 = lim_dict[vn][1]
    ax.plot([v0, v1], [v0, v1], '-k')
    ax.text(.05, .9, '('+abc[pp-1]+') '+vn_dict[vn], fontweight='bold',
        transform=ax.transAxes, size=fs, bbox=dict(facecolor='w', edgecolor='None', alpha=balpha))
    ax.text(.05, .8, errstat(Bc, vn), transform=ax.transAxes, size=fs, va='top',
        bbox=dict(facecolor='w', edgecolor='None', alpha=balpha))
    if pp in [3,4]:
        ax.set_xlabel('Observations', size=fs)
    else:
        ax.set_xlabel('')
    if pp in [1,3]:
        ax.set_ylabel('Model', size=fs)
    else:
        ax.set_ylabel('')
    if pp==1:
        ax.set_xticks([15,20,25,30,35])
        ax.set_yticks([15,20,25,30,35])
        ax.text(.95, .24, 'Z =    0 m', color='r', fontweight='bold',
            transform=ax.transAxes, horizontalalignment='right', size=fs)
        ax.text(.95, .17, 'Z = -10 m', color='g', fontweight='bold',
            transform=ax.transAxes, horizontalalignment='right', size=fs)
        ax.text(.95, .1, 'Z = -30 m', color='b', fontweight='bold',
            transform=ax.transAxes, horizontalalignment='right', size=fs)
        ax.text(.95, .03, 'Z = -80 m', color='k', fontweight='bold',
            transform=ax.transAxes, horizontalalignment='right', size=fs)
            
    ax.tick_params(labelsize=.8*fs)
    
    pp+=1

fig.tight_layout()
plt.savefig(out_plot_fn)
plt.show()