"""
Plot cast data for WOAC cruises, comparing with a LiveOcean run
"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun

# choose whether to properties vs. z or property-property scatter plots
plot_type = 'scatter' # 'casts' or 'scatter'
testing = False

# here is where we specify the ROMS run to use
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
Ldir['gtagex'] = Ldir['gtag']+'_'+'lo8b'
year = 2017

import zfun

pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np

import netCDF4 as nc

# load data
in_dir = '../../ptools_data/woac/'
Casts = pd.read_pickle(in_dir + 'Casts_'+str(year)+'.p')
sta_df = pd.read_pickle(in_dir + 'sta_df.p')

# prepare for output
out_dir0 = '../../ptools_output/woac/'
Lfun.make_dir(out_dir0)
out_dir = out_dir0 + 'val_' + plot_type + '_' + Ldir['gtagex'] + '/'
Lfun.make_dir(out_dir, clean=True)

# PLOTTING
plt.close('all')

# variable names in Casts, for reference
# ['Temp. (deg C)', 'Salinity', 'Sigma (kg m-3)', 'DO Flag', 'NO3 (uM)',
#        'NO2 (uM)', 'NH4 (uM)', 'Chl (ug L-1)', 'TA1 (umol kg-1)',
#        'DIC1 (umol kg-1)', 'TA1 Flag', 'DIC1 Flag', 'TA2 (umol kg-1)',
#        'DIC2 (umol kg-1)', 'TA2 Flag', 'DIC2 Flag', 'pH', 'pCO2 (uatm)',
#        'CO2 (umol kg-1)', 'HCO3- (umol kg-1)', 'CO3-- (umol kg-1)', 'Omega Ca',
#        'Omega Ar', 'DO (uM)']
              
# clean up data with bad flags
Casts.loc[Casts['DO Flag']==3,'DO (uM)'] = np.nan
Casts.loc[Casts['DO Flag']==4,'DO (uM)'] = np.nan
Casts.loc[Casts['TA1 Flag']==3,'TA1 (umol kg-1)'] = np.nan
Casts.loc[Casts['TA1 Flag']==4,'TA1 (umol kg-1)'] = np.nan
Casts.loc[Casts['DIC1 Flag']==3,'DIC1 (umol kg-1)'] = np.nan
Casts.loc[Casts['DIC1 Flag']==4,'DIC1 (umol kg-1)'] = np.nan

# make some derived variables
Casts['DIN (uM)'] = Casts['NO3 (uM)'] + Casts['NO2 (uM)'] + Casts['NH4 (uM)']
Casts['DIC (uM)'] = Casts['DIC1 (umol kg-1)'] * (1000 + Casts['Sigma (kg m-3)'])/1000
Casts['TA (uM)'] = Casts['TA1 (umol kg-1)'] * (1000 + Casts['Sigma (kg m-3)'])/1000

# list of variables to plot
vn_list = ['Temp. (deg C)', 'Salinity', 'DIN (uM)', 'DO (uM)', 'TA (uM)', 'DIC (uM)']
# associated model names
#vn_list_mod = ['temp', 'salt', 'NO3', 'oxygen', 'TIC', 'alkalinity']
vn_list_mod = ['temp', 'salt', 'NO3', 'oxygen', 'alkalinity', 'TIC']
vn_dict = dict(zip(vn_list, vn_list_mod))
# and plot panel numbers
np_dict = dict(zip(vn_list,[1, 2, 5, 6, 9, 10]))
       
# limits for plotting
lim_dict = {'temp': (0, 20),
        'salt': (15, 35),
        'NO3': (0, 48),
        'phytoplankton': (-1, 14),
        'oxygen': (0, 450),
        'TIC': (1000, 2600),
        'alkalinity': (1000, 2600)}

plt.close('all')

if testing == True:
    cn_list = [0]
else:
    cn_list = sta_df.index
    
for cn in cn_list:
    
    # get the corresponding model cast
    cinfo = sta_df.loc[cn,:]
    in_dir_m = Ldir['LOo'] + 'cast/'+Ldir['gtagex']+'/'
    fn_m = (in_dir_m + 'WOAC' + str(cinfo['Station']) + '_' +
            cinfo['Datetime'].strftime('%Y.%m.%d') + '.nc')
    ds = nc.Dataset(fn_m)
    
    fig = plt.figure(figsize=(22,11))
    
    c = Casts[Casts['castnum']==cn]
    
    # drop duplicate depths, and sort by Z
    # (only needed for a few RBTSON201705 casts)
    cc = c.set_index('Z (m)')
    cc['Z (m)'] = cc.index
    cc = cc[~cc.index.duplicated()]
    cc = cc.sort_index() # sort deepest to shallowest
    
    if plot_type == 'casts':    
        for vn in vn_list:
            ax = fig.add_subplot(3,4,np_dict[vn])
            try:
                cc.plot(x=vn, y='Z (m)', ax=ax, legend=False, style='-*b')
            except ValueError:
                pass
            ax.text(.05,.1, vn, transform=ax.transAxes, fontweight='bold')
            ax.set_xlabel('')
            # add model line
            mod_vn = vn_dict[vn]
            ax.plot(ds[mod_vn][:], ds['z_rho'], '-r')
            # set limits
            ax.set_xlim(lim_dict[mod_vn])
            ax.set_ylim(-200,0)
            ax.grid(True)
        
            if vn == 'Temp. (deg C)':
                ax.text(.05, .9, 'Observed', color='b', transform=ax.transAxes, fontweight='bold')
                ax.text(.05, .8, 'Modeled', color='r', transform=ax.transAxes, fontweight='bold')
                
    elif plot_type == 'scatter':
        
        # vn_list = ['Temp. (deg C)', 'Salinity', 'DIN (uM)', 'DO (uM)', 'TA (uM)', 'DIC (uM)']

        vn_pairs = [('Salinity', 'Temp. (deg C)'), ('TA (uM)', 'DIC (uM)'), ('DO (uM)', 'DIN (uM)'), ('Salinity', 'TA (uM)')]
        # and plot panel numbers
        np_pair_dict = dict(zip(vn_pairs,[1, 2, 5, 6]))
        
        for vn_pair in vn_pairs:
        
            ax = fig.add_subplot(2,4,np_pair_dict[vn_pair])
            try:
                x = cc[vn_pair[0]].values
                y = cc[vn_pair[1]].values
                z = cc['Z (m)'].values
                
                mod_vnx = vn_dict[vn_pair[0]]
                mod_vny = vn_dict[vn_pair[1]]
                
                modx = ds[mod_vnx][:]
                mody = ds[mod_vny][:]
                modz = ds['z_rho'][:]
                
                zz_list = [(-300,-60), (-60, -30), (-30,-10), (-10,0)]
                zc_dict = dict(zip(zz_list, ['darkblue', 'green', 'gold', 'red']))
                for zz in zz_list:
                    z0 = zz[0]
                    z1 = zz[1]
                    mask = (z > z0) & (z <= z1)
                    ax.plot(x[mask], y[mask], marker='*', color=zc_dict[zz], markersize=14, linestyle='')
                    
                    modmask = (modz > z0) & (modz <= z1)
                    ax.plot(modx[modmask], mody[modmask], marker='o', color=zc_dict[zz], markersize=5, alpha=.5, linestyle='')
                    
                    
                if vn_pair == ('Salinity', 'Temp. (deg C)'):
                    counter = 0
                    for zz in zz_list:
                        z0 = zz[0]
                        z1 = zz[1]
                        ax.text(.05, .6+.06*counter, 'Depth '+str(-z1)+'-'+str(-z0)+' m', color=zc_dict[zz], transform=ax.transAxes, fontweight='bold')
                        counter += 1
                    ax.text(.05, .2, 'Star = Observed', color='k', transform=ax.transAxes, fontweight='bold')
                    ax.text(.05, .1, 'Circle = Modeled', color='k', transform=ax.transAxes, fontweight='bold', alpha=.5)
                        
            except ValueError:
                print('missing data for ' + str(vn_pair)) 
                pass
            #ax.text(.05,.1, vn, transform=ax.transAxes, fontweight='bold')
            ax.set_xlabel(vn_pair[0])
            ax.set_ylabel(vn_pair[1])
            # set limits
            ax.set_xlim(lim_dict[mod_vnx])
            ax.set_ylim(lim_dict[mod_vny])
            ax.grid(True)
    
    if testing == True:
        pass
    else:
        ds.close()
    
    # add location map
    ax = fig.add_subplot(1,2,2)
    pfun.add_coast(ax)
    ax.axis([-125.5, -122, 47, 49])
    pfun.dar(ax)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    
    cr = cc['Cruise'].values[0]
    sta_df_cr = sta_df[sta_df['Cruise'] == cr]
    
    for cnn in sta_df_cr.index:
        x = sta_df_cr.loc[cnn, 'Longitude']
        y = sta_df_cr.loc[cnn, 'Latitude']
        ss = str(sta_df_cr.loc[cnn, 'Station'])
        ax.plot(x,y, '*b', alpha=.3)
        ax.text(x, y+.02, ss, fontsize=15, alpha = .3)
        
    x = sta_df.loc[cn, 'Longitude']
    y = sta_df.loc[cn, 'Latitude']
    ss = str(sta_df.loc[cn, 'Station'])
    ax.plot(x,y, '*b', markersize=20)
    ax.text(x, y+.02, ss, fontsize=15)
    
    ttext = ('Cruise = %s, Station = %s, Date = %s' % (sta_df.loc[cn,'Cruise'],
        sta_df.loc[cn,'Station'],
        datetime.strftime(sta_df.loc[cn,'Datetime'],'%Y.%m.%d')))
    
    ax.set_title(ttext)
    
    if testing == True:
        plt.show()
    else:
        plt.savefig(out_dir + sta_df.loc[cn,'Cruise'] + '_' +
            str(sta_df.loc[cn,'Station']) + '_' +
            datetime.strftime(sta_df.loc[cn,'Datetime'],'%Y.%m.%d') + '.png')
        plt.close()
        
       
    
    