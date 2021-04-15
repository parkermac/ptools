"""
Plot WOAC stations in an informative way.

"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3', ex_name='lo8b')
import plotting_functions as pfun

# get the list of cruises and casts
sta_df_dir = Path(__file__).absolute().parent.parent.parent / 'ptools_output' / 'woac2'
OMa = pd.read_pickle(sta_df_dir / 'OMa.p')
sta_df = pd.read_pickle(sta_df_dir / 'sta_df.p')

cruises = sta_df['Cruise'].unique()

fs=16
plt.rc('font', size=fs)
plt.close('all')
fig = plt.figure(figsize=(14,11))
ii = 1

# sort cruises by date
cruise_list = cruise_list = ['CAB1079', 'AQ201710', 'RC001',
    'RBTSN201805', 'RC006', 'RC007', 'NORSEMANIIOCT2018']

for cruise in cruise_list:
    ax = fig.add_subplot(3,3,ii)
    
    c_df = sta_df[sta_df['Cruise']==cruise]
    
    c_df.plot(x='Longitude', y = 'Latitude', style='or', alpha=.5,
        ax=ax, legend=False)
    for sn in c_df.index:
        ax.text(c_df.loc[sn,'Longitude']+.02, c_df.loc[sn,'Latitude'],str(int(sn)), size=10)
    pfun.add_coast(ax, color='gray')
    pfun.dar(ax)
    ax.axis([-125, -122, 47, 48.5])
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks([-124, -123])
    ax.set_yticks([47, 48])
    ax.set_title(cruise)
    
    dt = c_df['Datetime'].mean()
    ax.text(.1,.1,'%d/%d' % (dt.month, dt.year), transform=ax.transAxes, weight='bold')
    
    
    ii += 1
    
fig.tight_layout()
plt.show()
plt.rcdefaults()