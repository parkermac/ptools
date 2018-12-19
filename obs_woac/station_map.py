"""
Plot a station map
"""
import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()

pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

import pandas as pd
import matplotlib.pyplot as plt

# load data
in_dir = '../../ptools_data/woac/'
Casts = pd.read_pickle(in_dir + 'Casts_2017.p')
sta_df = pd.read_pickle(in_dir + 'sta_df.p')

# prepare for output
out_dir = '../../ptools_output/woac/'
Lfun.make_dir(out_dir)

# PLOTTING
plt.close('all')

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)

# This plots ALL the data locations, to make sure that our averaging into
# single lat,lon values for a given station is appropriate.
Casts.plot(x='Longitude', y='Latitude',
    style='*c', alpha=.3, ax = ax, legend=False)
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis([-125.5, -122, 47, 49])
ax.set_title('WOAC Stations 2017')
ax.set_xlabel('Longitude (deg)')
ax.set_ylabel('Latitude (deg)')

cruise_list = list(set(sta_df['Cruise']))
color_list = ['r','g','b']
color_list = color_list[:len(cruise_list)]
color_dict = dict(zip(cruise_list, color_list))
offset_list = [0,.02,-.02]
offset_list = offset_list[:len(cruise_list)]
offset_dict = dict(zip(cruise_list, offset_list))

for cn in sta_df.index:
    x = sta_df.loc[cn, 'Longitude']
    y = sta_df.loc[cn, 'Latitude']
    ss = str(sta_df.loc[cn, 'Station'])
    cr = sta_df.loc[cn,'Cruise']
    #print('%5s: x = %8.2f y = %8.2f' % (ss,x,y))
    ax.text(x, y+offset_dict[cr], ss, fontsize=15, color=color_dict[cr], fontweight='bold')
    
ii = 0
for cr in color_dict.keys():
    ax.text(.05, .1+ii*.1, cr, fontsize=20, color=color_dict[cr], transform=ax.transAxes)
    ii+=1

plt.show()

plt.savefig(out_dir + 'station_map.png')