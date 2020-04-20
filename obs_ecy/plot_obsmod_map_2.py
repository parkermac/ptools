"""
Plots data from the Department of Ecology, combining data from CTD casts
and bottles at the same station.  Also including model output.

This makes a summary map, using color in side-by-side boxes.

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

import matplotlib
import matplotlib.cm as cm

v_dict = {
    'Salinity': (20,34,'Salinity','salt',[0, -10, -80]),
    'Temp. (deg C)': (8,18,'Temp. [$^{\circ} C$]','temp',[0, -10, -80]),
    'DO (mg L-1)': (0,12,'DO [$mg\ L^{-1}$]','DO',[0, -10, -80]),
    'DIN (uM)': (0,35,'DIN [$\mu M$]','DIN',[0, -10, -30])
    }
    
v_list = list(v_dict.keys())
if testing:
    v_list = ['Salinity']

# PLOTTING
plt.close('all')
fs = 18
dx = .04
sz = 10
abc = 'abcd'
for vname in v_list:

    fig = plt.figure(figsize=(20,12))
    v0, v1, vname_str, short_name, z_list = v_dict[vname]
    
    norm = matplotlib.colors.Normalize(vmin=v0, vmax=v1, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cm.rainbow)
    
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
            
            if np.isnan(vo):
                voc='None'
            else:
                voc = mapper.to_rgba(vo)

            if np.isnan(vm):
                vmc='None'
            else:
                vmc = mapper.to_rgba(vm)
    
            ax.plot(xs-dx, ys, 's', markerfacecolor=voc, markeredgecolor='k',
                markersize=sz)
            ax.plot(xs+dx, ys, 's', markerfacecolor=vmc, markeredgecolor='k',
                markersize=sz)
            
            # if counter == 3:
            #     ax.text(xs, ys, sn, ha='center', va='center', size=.5*fs)
        
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
            lx = -122.5; ly = 46.7; ldx = .2; lsz = 50
            ax.plot(lx-ldx, ly, 's', markersize=lsz, markerfacecolor='None', markeredgecolor='k')
            ax.text(lx-ldx, ly, 'Obs', size=fs, ha='center', va='center')
            ax.plot(lx+ldx, ly, 's', markersize=lsz, markerfacecolor='None', markeredgecolor='k')
            ax.text(lx+ldx, ly, 'Mod', size=fs, ha='center', va='center')
            
            # Inset colorbar
            from mpl_toolkits.axes_grid1.inset_locator import inset_axes
            cbaxes = inset_axes(ax, width="4%", height="40%", loc=2, borderpad=2.5) 
            cb = fig.colorbar(mapper, cax=cbaxes, orientation='vertical')
            
            plt.setp(cb.ax.yaxis.get_ticklabels(), fontsize=fs)
            
        counter += 1
    
    outname = Ldir['gtagex'] + '_' + year_str + '_' + short_name + '.png'
    out_fn = outdir + outname
    #fig.tight_layout()
    if not testing:
        plt.savefig(out_fn)
plt.show()



