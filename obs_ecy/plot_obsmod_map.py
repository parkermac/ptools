"""
Plots data from the Department of Ecology, combining data from CTD casts
and bottles at the same station.  Also including model output.

This makes a summary map.  Now the AREA of symbols is proportional to
the value they represent (had been the radius).

"""
import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
import zfun

sys.path.append(os.path.abspath('../../LiveOcean/plotting'))
import pfun

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seawater as sw

testing = False

Ldir = Lfun.Lstart()
Ldir['gtagex'] = 'cas6_v3_lo8b'
year = 2017
year_str = str(year)

# input
dir0 = Ldir['parent'] + 'ptools_output/ecology/'
out_fn = 'ObsMod_' + Ldir['gtagex'] + '_'+ str(year) + '.p'
Bc = pd.read_pickle(dir0 + out_fn)

# output
outdir = dir0 + 'obsmod_maps/'
Lfun.make_dir(outdir)

# Bc is a DataFrame with processed bottles and casts, both observed and modeled
# at specific depths:
#      Station       Date  Znom    Z  ... Mod DO (mg L-1) Mod DIN (uM) Density (kg m-3) Mod Density (kg m-3)
# 0     ADM001 2017-02-23     0  NaN  ...         10.0564      16.9518          1022.36              1021.47
# 1     ADM001 2017-02-23   -10  NaN  ...          9.2235      20.4011          1022.59              1022.89
# 2     ADM001 2017-02-23   -30  NaN  ...         8.91924      21.2567          1022.67              1023.31
# 3     ADM001 2017-02-23   -80  NaN  ...         8.79093      20.8736          1023.34              1023.72
# The full ist of columns is:
# Index(['Station', 'Date', 'Znom', 'Z', 'Salinity', 'Temp. (deg C)',
#        'Chl (mg m-3)', 'DO (mg L-1)', 'DIN (uM)', 'Mod Salinity',
#        'Mod Temp. (deg C)', 'Mod Chl (mg m-3)', 'Mod DO (mg L-1)',
#        'Mod DIN (uM)', 'Density (kg m-3)', 'Mod Density (kg m-3)'],
#       dtype='object')

Bc.index = Bc.Date

dir_ecy = '../../ptools_data/ecology/'
# load processed station info and data
sta_df = pd.read_pickle(dir_ecy + 'sta_df.p')
# add Canadian data
dir_ca = '../../ptools_data/canada/'
# load processed station info and data
sta_df_ca = pd.read_pickle(dir_ca + 'sta_df.p')
# merge
sta_df = pd.concat((sta_df, sta_df_ca), sort=False)

# axis limits
x0 = -124.5; x1 = -122; y0 = 46; y1 = 49.5 # Salish Sea
aa = [x0, x1, y0, y1]

# def msz(vv, v0, v1, scl=40):
#     msz = scl*(vv - v0)/(v1 - v0)
#     if msz < 1:
#         msz = 6.3/40
#     msz_new = (40/6.3)*np.sqrt(msz)
#     return msz_new
def msz(vv, v0, v1):
    vvs = (vv - v0)/(v1 - v0)
    scl = 40
    if vvs < 0:
        vvs = 1/scl**2
    msz = scl*np.sqrt(vvs)
    return msz

v_dict = {
    'Salinity': (16,34,'Salinity','salt',[0, -10, -80]),
    'Temp. (deg C)': (8,18,'Temp. [$^{\circ} C$]','temp',[0, -10, -80]),
    'DO (mg L-1)': (3,10,'DO [$mg\ L^{-1}$]','DO',[0, -10, -80]),
    'DIN (uM)': (0,35,'DIN [$\mu M$]','DIN',[0, -10, -30])
    }
    
v_list = list(v_dict.keys())
if testing:
    v_list = ['Salinity']

# PLOTTING
plt.close('all')
fs = 18
co = 'c'; ao = .5
cm = 'b'
abc = 'abcd'
for vname in v_list:

    fig = plt.figure(figsize=(20,12))
    v0, v1, vname_str, short_name, z_list = v_dict[vname]
    v01 = (v0+v1)/2
    
    counter = 1
    for z in z_list:
    
        ax = fig.add_subplot(1, 3, counter)
        pfun.add_coast(ax, color='gray')
        ax.axis(aa)
        pfun.dar(ax)
        ax.tick_params(labelsize=.8*fs)
        ax.set_xlabel('Longitude', size=fs)
        
        if counter > 1:
            ax.set_yticklabels([])

        for sn in sta_df.index:
            xs = sta_df.loc[sn, 'Longitude']
            ys = sta_df.loc[sn, 'Latitude']
    
            mo0 = 6; mo1 = 9 # month range to consider
            mask = ((Bc.Station==sn) & (Bc.Znom==z)
                & (Bc.index.month >= mo0) & (Bc.index.month <= mo1))
    
            vo = Bc.loc[mask,vname].mean()
            vm = Bc.loc[mask,'Mod ' + vname].mean()
            
            # print('Obs. Value = %0.1f msz = %0.1f' % (vo, msz(vo, v0, v1)))
            # print('Mod. Value = %0.1f msz = %0.1f' % (vm, msz(vm, v0, v1)))
    
            ax.plot(xs, ys, 'o', markerfacecolor=co, markeredgecolor=co,
                markersize=msz(vo, v0, v1), alpha=ao)
            ax.plot(xs, ys, 'o', markerfacecolor='None', markeredgecolor=cm,
                markersize=msz(vm, v0, v1))
            
            if counter == 3:
                ax.text(xs, ys, sn, ha='center', va='center', size=.5*fs)
        
        if counter == 1:
            ax.text(.95, .98, '(%s) Z = %d (m)\nMonths %d to %d' % (abc[counter-1], z, mo0, mo1),
                ha='right', va='top', size=fs, transform=ax.transAxes,
                bbox=dict(facecolor='w', edgecolor='None'))
            ax.text(.95, .9, vname_str, weight='bold',
                ha='right', va='top', size=1.3*fs, transform=ax.transAxes,
                bbox=dict(facecolor='w', edgecolor='None'))
            
        else:
            ax.text(.95, .98, '(%s) Z = %d (m)' % (abc[counter-1], z),
                ha='right', va='top', size=fs, transform=ax.transAxes,
                bbox=dict(facecolor='w', edgecolor='None'))
    
        if counter == 1:
            ax.set_ylabel('Latitude', size=fs)
            
            # legend
            lx = .7; ly = .15
            ax.plot(lx, ly, 'o', markerfacecolor=co, markeredgecolor=co,
                markersize=msz(v1, v0, v1), alpha=ao, transform=ax.transAxes)
            ax.plot(lx, ly, 'o', markerfacecolor='None', markeredgecolor=cm,
                markersize=msz(v0, v0, v1), transform=ax.transAxes)
            ax.plot(lx, ly, 'o', markerfacecolor='None', markeredgecolor=cm,
                markersize=msz(v01, v0, v1), transform=ax.transAxes)
            ax.plot(lx, ly, 'o', markerfacecolor='None', markeredgecolor=cm,
                markersize=msz(v1, v0, v1), transform=ax.transAxes)
            ax.text(lx, ly+.05, 'Observed (%d)' % (v1),
                ha='center', va='center', size=fs, style='italic', color=co,
                transform=ax.transAxes)
            ax.text(lx, ly-.06, 'Modeled (%d, %d, %d)' % (v0, v01, v1),
                ha='center', va='center', size=fs, style='italic', color=cm,
                transform=ax.transAxes)
            ax.text(lx, .03, '$Circle\ Area\ \propto\ Value$',
                ha='center', va='center', size=.8*fs, style='italic', color='k',
                transform=ax.transAxes,bbox=dict(facecolor='w', edgecolor='None'))
                        
        counter += 1
    
    outname = Ldir['gtagex'] + '_' + year_str + '_' + short_name + '_CIRCLES.png'
    out_fn = outdir + outname
    #fig.tight_layout()
    if not testing:
        plt.savefig(out_fn)
plt.show()



