"""
This plots results of obs_mod_compare.py as property=property "pp" plots.

"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

pth = Path(__file__).absolute().parent.parent.parent / 'LiveOcean' / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun
Ldir = Lfun.Lstart()
out_dir = Ldir['parent'] + 'ptools_output/woac2_plots/'
Lfun.make_dir(out_dir)

pth = Path(__file__).absolute().parent.parent.parent / 'LiveOcean' / 'plotting'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import pfun

# get the list of cruises and casts
sta_df_dir = Path(__file__).absolute().parent.parent.parent / 'ptools_output' / 'woac2'
OMa = pd.read_pickle(sta_df_dir / 'OMa.p')
sta_df = pd.read_pickle(sta_df_dir / 'sta_df.p')
sta_df['Station'] = sta_df.index

# Names of variables
# vv_list = ['s', 'th', 'do', 'din', 'ta', 'dic', 'ph', 'arag']
vn_text = {'s': 'Salinity',
            'th': 'Potential Temperature (degC)',
        'do': 'Dissolved Oxygen (umol/kg)',
        'din': 'Dissolved Inorganic Nitrogen (umol/kg)',
        'ta': 'Alkalinity (umol/kg)',
        'dic': 'Dissolved Inorganic Carbon (umol/kg)',
        'ph': 'pH',
        'arag': 'Aragonite Saturation State',
    }


# cruises mostly in Puget Sound
cruise_list_ps = cruise_list = ['CAB1079', 'RC001', 'RC006', 'RC007']
    
# cruises mostly in Juan de Fuca
cruise_list_jdf = cruise_list = ['AQ201710', 'RBTSN201805', 'NORSEMANIIOCT2018']

# list of stations in Hood Canal
hcs = [11, 12, 13, 14, 15, 401, 402]


fs=14
plt.rc('font', size=fs)
plt.close('all')
alpha=.5

for vn in ['s', 'th', 'do', 'din', 'ta', 'dic', 'ph', 'arag']:

    fig = plt.figure(figsize=(16,10))
    ii = 3
    for cruise in cruise_list_ps:
        c_df = sta_df[sta_df['Cruise']==cruise]
    
        if cruise == cruise_list_ps[1]:
            ax=fig.add_subplot(2,3,1)
            c_df[c_df['Station'].isin(hcs)].plot(x='Longitude', y = 'Latitude', style='*k', alpha=alpha, markersize=12,
                ax=ax, legend=False)
            c_df[~c_df['Station'].isin(hcs)].plot(x='Longitude', y = 'Latitude', style='ok', alpha=alpha, markersize=6,
                ax=ax, legend=False)
            pfun.add_coast(ax, color='gray')
            pfun.dar(ax)
            ax.axis([-125, -122, 47, 48.5])
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_xticks([-124, -123])
            ax.set_yticks([47, 48])
    
        om = OMa[OMa['Cruise']==cruise]
        vnm = vn + '_m'
        ax = fig.add_subplot(2,3,ii)
        vmin = np.min((OMa[vn].min(),OMa[vnm].min()))
        vmax = np.max((OMa[vn].max(),OMa[vnm].max()))
    
        om1 = om[om['Station'].isin(hcs)]
        om2 = om[~om['Station'].isin(hcs)]
    
        om1[om1['z']>=-20].plot(x=vn, y=vnm, legend=False, xlim=(vmin,vmax), ylim=(vmin,vmax),
            ax=ax, style='*r',markersize=12, alpha=alpha)
        om1[om1['z']<-20].plot(x=vn, y=vnm, legend=False, xlim=(vmin,vmax), ylim=(vmin,vmax),
            ax=ax, style='*b',markersize=12, alpha=alpha)

        om2[om2['z']>=-20].plot(x=vn, y=vnm, legend=False, xlim=(vmin,vmax), ylim=(vmin,vmax),
            ax=ax, style='.r',markersize=6, alpha=alpha)
        om2[om2['z']<-20].plot(x=vn, y=vnm, legend=False, xlim=(vmin,vmax), ylim=(vmin,vmax),
            ax=ax, style='.b',markersize=6, alpha=alpha)
        
        dt = c_df['Datetime'].mean()
        ax.text(.95,.05,'%d/%d' % (dt.month, dt.year), transform=ax.transAxes, weight='bold', ha='right')
        ax.text(.05,.95,cruise, transform=ax.transAxes, weight='bold', va='top')

        ax.plot([vmin,vmax],[vmin,vmax],'-k')
        ax.axis('square')
        ax.set_xlabel('')
        ax.set_ylabel('')
    
        ii += 1
    fig.suptitle('PS Stations (x=obs, y=mod) RED = above 20 m: ' + vn_text[vn])
    fig.savefig(out_dir + vn + '_ps.png')
    
    fig = plt.figure(figsize=(12,10))
    ii = 2
    for cruise in cruise_list_jdf:
        c_df = sta_df[sta_df['Cruise']==cruise]
    
        if cruise == cruise_list_jdf[1]:
            ax=fig.add_subplot(2,2,1)
            c_df.plot(x='Longitude', y = 'Latitude', style='ok', alpha=alpha, markersize=6,
                ax=ax, legend=False)
            pfun.add_coast(ax, color='gray')
            pfun.dar(ax)
            ax.axis([-125, -122, 47, 48.5])
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_xticks([-124, -123])
            ax.set_yticks([47, 48])
    
        om = OMa[OMa['Cruise']==cruise]
        vnm = vn + '_m'
        ax = fig.add_subplot(2,2,ii)
        vmin = np.min((OMa[vn].min(),OMa[vnm].min()))
        vmax = np.max((OMa[vn].max(),OMa[vnm].max()))
    
        om[om['z']>=-20].plot(x=vn, y=vnm, legend=False, xlim=(vmin,vmax), ylim=(vmin,vmax),
            ax=ax, style='.r',markersize=6, alpha=alpha)
        om[om['z']<-20].plot(x=vn, y=vnm, legend=False, xlim=(vmin,vmax), ylim=(vmin,vmax),
            ax=ax, style='.b',markersize=6, alpha=alpha)
        
        dt = c_df['Datetime'].mean()
        ax.text(.95,.05,'%d/%d' % (dt.month, dt.year), transform=ax.transAxes, weight='bold', ha='right')
        ax.text(.05,.95,cruise, transform=ax.transAxes, weight='bold', va='top')

        ax.plot([vmin,vmax],[vmin,vmax],'-k')
        ax.axis('square')
        ax.set_xlabel('')
        ax.set_ylabel('')
    
        ii += 1
    fig.suptitle('JdF Stations (x=obs, y=mod) RED = above 20 m: ' + vn_text[vn])
    fig.savefig(out_dir + vn + '_jdf.png')
    
#plt.show()
plt.rcdefaults()
