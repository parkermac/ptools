"""
Plot cast data.
"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()
import zfun

pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# load data
in_dir = '../../ptools_data/woac/'
Casts = pd.read_pickle(in_dir + 'Casts_2017.p')
sta_df = pd.read_pickle(in_dir + 'sta_df.p')

# prepare for output
out_dir = '../../ptools_output/woac/'
Lfun.make_dir(out_dir)

# PLOTTING
plt.close('all')

vn_list = ['Temp. (deg C)', 'Salinity', 'Sigma (kg m-3)', 'DO Flag', 'NO3 (uM)',
       'NO2 (uM)', 'NH4 (uM)', 'Chl (ug L-1)', 'TA1 (umol kg-1)',
       'DIC1 (umol kg-1)', 'TA1 Flag', 'DIC1 Flag', 'TA2 (umol kg-1)',
       'DIC2 (umol kg-1)', 'TA2 Flag', 'DIC2 Flag', 'pH', 'pCO2 (uatm)',
       'CO2 (umol kg-1)', 'HCO3- (umol kg-1)', 'CO3-- (umol kg-1)', 'Omega Ca',
       'Omega Ar', 'DO (uM)']

NR, NC = zfun.get_rc(len(vn_list))

for cn in sta_df.index[:10]:
    
    fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(14,10), squeeze=False)
    
    c = Casts[Casts['castnum']==cn]
    
    ii = 0
    for vn in vn_list:
        ir, ic = zfun.get_irc(ii, NC)
        ax = axes[ir,ic]
        c.plot(x=vn, y='Z (m)', ax=ax, legend=False, style='-*b')
        ax.text(.05,.1, vn, transform=ax.transAxes, fontweight='bold')
        
        ii += 1
    
    ttext = ('Cruise = %s, Station = %s, Date = %s' % (sta_df.loc[cn,'Cruise'],
        sta_df.loc[cn,'Station'],
        datetime.strftime(sta_df.loc[cn,'Datetime'],'%Y.%m.%d')))
    fig.suptitle(ttext)
        
plt.show()
       
    
    