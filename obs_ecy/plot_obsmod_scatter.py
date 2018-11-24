#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots data from the Department of Ecology, combining data from CTD casts
and bottles at the same station.  Also including model output.

This makes the summary scatterplot for all variables of interest.

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

Ldir['gtagex'] = 'cas4_v2_lo6biom'
year = 2017
dir11 = Ldir['parent'] + 'ptools_output/ecology/'
out_fn = 'ObsMod_' + Ldir['gtagex'] + '_'+ str(year) + '.p'
print('Loading ' + out_fn)
Bc = pd.read_pickle(dir11 + out_fn)

# add density to calculate stratification
Bc['Density (kg m-3)'] = sw.dens0(Bc['Salinity'], Bc['Temp. (deg C)'])
Bc['Mod Density (kg m-3)'] = sw.dens0(Bc['Mod Salinity'], Bc['Mod Temp. (deg C)'])

D0 = Bc[Bc['Znom']==0]
D10 = Bc[Bc['Znom']==-10]
D30 = Bc[Bc['Znom']==-30]

S5 = D10['Density (kg m-3)'].values - D0['Density (kg m-3)'].values
S20 = D30['Density (kg m-3)'].values - D10['Density (kg m-3)'].values
SM5 = D10['Mod Density (kg m-3)'].values - D0['Mod Density (kg m-3)'].values
SM20 = D30['Mod Density (kg m-3)'].values - D10['Mod Density (kg m-3)'].values


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
    if pp in [4,5,6]:
        ax.set_xlabel('Observations')
    else:
        ax.set_xlabel('')
    if pp in [1,4]:
        ax.set_ylabel('Model')
    else:
        ax.set_ylabel('')
    if pp==1:
        ax.text(.95, .24, 'Z =   0 m', color='r', fontweight='bold',
            transform=ax.transAxes, horizontalalignment='right')
        ax.text(.95, .17, 'Z = -10 m', color='g', fontweight='bold',
            transform=ax.transAxes, horizontalalignment='right')
        ax.text(.95, .1, 'Z = -30 m', color='b', fontweight='bold',
            transform=ax.transAxes, horizontalalignment='right')
    
    pp+=1

# generate error stats
SS = pd.DataFrame(columns=['Strat. (kg m-3)', 'Mod Strat. (kg m-3)'])
SS.loc[:,'Strat. (kg m-3)'] = np.concatenate((S5,S20))
SS.loc[:,'Mod Strat. (kg m-3)'] = np.concatenate((SM5,SM20))

# plot the stratification
ax = fig.add_subplot(NR, NC ,6)
ax.plot(S5, SM5, 'om', alpha=.4)
ax.plot(S20, SM20, 'oc', alpha=.4)
v0 = 0; v1 = 3
ax.axis([v0, v1, v0, v1])
ax.grid(True)
ax.plot([v0, v1], [v0, v1], '-k')
ax.set_title(errstat(SS, 'Strat. (kg m-3)'))
ax.set_xlabel('Observations')
#ax.set_ylabel('Model')
ax.text(.05, .9, 'Stratification (kg m-3)', fontweight='bold', transform=ax.transAxes)
ax.text(.95, .24, '0-10 m', color='m', fontweight='bold',
    transform=ax.transAxes, horizontalalignment='right')
ax.text(.95, .17, '10-30 m', color='c', fontweight='bold',
    transform=ax.transAxes, horizontalalignment='right')

fig.tight_layout()
plt.show()