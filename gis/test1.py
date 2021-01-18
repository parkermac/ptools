"""
First code to try reading some GIS files from Ecolgy using geopandas.
"""

import os, sys

sys.path.append(os.path.abspath('../../LiveOcean/plotting'))
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import pfun

import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd

s_fn = '../../ptools_data/Ecology_Point_Sources_2019_10/SSM2_inflows_point.shp'

s = gpd.read_file(s_fn)
# s.crs gives projection info
# convert to lon, lat:
ss = s.to_crs('EPSG:4326') # apparently 4326 is the magic number!

# we can also read in time series for each station
id = '566'
sn = 'Bellingham'
w_fn = '../../ptools_data/Ecology_Point_Sources_2019_10/WWTPzipfiles/'+id+'_'+sn+'.xlsx'
w = pd.read_excel(w_fn, skiprows=1)
w = w.set_index('Date')
# and this has these columns:
# Index(['Year', 'Month', 'Day', 'Hour', 'Minute', 'Bin1', 'Flow, cms',
#        'Temp (C)', 'Saln (ppt)', 'NH4 (mg/L)', 'NO3+NO2 (mg/L)', 'PO4 (mg/L)',
#        'DO (mg/L)', 'pH (SU)', 'DON(mg/L)', 'PON(mg/L)', 'DOP(mg/L)',
#        'POP(mg/L)', 'POCS(mg/L)', 'POCF(mg/L)', 'POCR(mg/L)', 'sDOC(mg/L)',
#        'fDOC (mg/L)', 'Diatoms', 'Dinoflag', 'Chlorophyll', 'DIC(mmoles/m3)',
#        'Alk(mmoles/m3)'],
#       dtype='object')

plt.close('all')
fig = plt.figure(figsize=(12,12))

ax = fig.add_subplot(111)
ss.plot(ax=ax, column='Depth', legend=True, cmap='rainbow')

pfun.add_coast(ax)

a = ss.total_bounds
pad = .1
aa = [a[0]-pad, a[2]+pad, a[1]-pad, a[3]+pad]
ax.axis(aa)

plt.show()