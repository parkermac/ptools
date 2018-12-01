"""
Process WOAC cruise data, and
save the results for future use.

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
import numpy as np
import pickle

import matplotlib.pyplot as plt

# where the data is, and where results will be stored
dir0 = '../../ptools_data/woac/'

a = pd.read_excel(dir0 + 'WOAC_data_9-7-2018_dataToParker.xlsx',
    parse_dates = ['Date_collected'])
    
# keep only selected columns
a = a[['record no', 'CRUISE_ID', 'Date_collected', 'Time_collected',
       'LONGITUDE_DEC', 'LATITUDE_DEC', 'STATION_NO', 'NISKIN_NO',
       'CTDPRS_DBAR', 'CTDTMP_DEG_C_ITS90', 'CTDSAL_PSS78', 'SIGMATHETA_KG_M3',
       'CTDOXY_UMOL_KG_ADJ', 'OXYGEN_FLAG_W', 'CTD pH ',
       'NITRATE_UMOL_L', 'NITRITE_UMOL_L', 'AMMONIUM_UMOL_L',
       'PHOSPHATE_UMOL_L', 'SILICATE_UMOL_L', 'CTD FLU (mg/m3)', 'CHLA (ug/l)',
       'CHLA 2 (ug/l)', 'CHLA avg (ug/l)', 'PHAEOPIGMENT (ug/l)',
       'PHAEOPIGMENT 2 (ug/l)', 'PHAEOPIGMENT avg (ug/l)', 'TA_1_umol_kg',
       'DIC_1_umol_kg', 'TA_1_QC ', 'DIC_1_QC', 'TA_2_umol_kg',
       'DIC_2_umol_kg', 'TA_2_QC ', 'DIC_2_QC', 'NITRATE umol_kg',
       'NITRITE umol_kg', 'AMMONIA umol_kg', 'PHOSPHATE umol_kg',
       'SILICATE umol_kg', 'Ph Total in situ',
       'pCO2 µatm', 'CO2 µmol/kg', 'HCO3- µmol/kg', 'CO3-- µmol/kg',
       'Omega Ca', 'Omega Ar']]

# get a list of station numbers
s = a.loc[:,'STATION_NO']
ss = set(s)
sl = []
for sta in ss:
    if isinstance(sta, int):
        sl.append(sta)
sl.sort()

# get a list of dates
d = a.loc[:,'Date_collected']
dd = list(set(d))
dd.sort()

# then associate station numbers with lat, lon
# go through casts
sta_df = pd.DataFrame(columns = ['Latitude', 'Longitude'])
sta_df.index.name = 'Station'
for date in dd:
    aa = a[a['Date_collected'] == date]
    for sta in sl:
        aaa = aa[aa['STATION_NO'] == sta]
        if len(aaa) > 0:
            lon = aaa.loc[:,'LONGITUDE_DEC'].mean()
            lat = aaa.loc[:,'LATITUDE_DEC'].mean()
            sta_df.loc[sta,'Latitude'] = lat
            sta_df.loc[sta,'Longitude'] = lon
            
# make a list of casts
for date in dd:
    aa = a[a['Date_collected'] == date]
    for sta in sl:
        aaa = aa[aa['STATION_NO'] == sta]
        if len(aaa) > 0:
            print('Date = ' + str(date) + ': Station = ' + str(sta))
    
# PLOTTING
plt.close('all')

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)

a.plot(x='LONGITUDE_DEC', y='LATITUDE_DEC', style='*r', alpha=.3, ax = ax)
pfun.dar(ax)
pfun.add_coast(ax)
ax.axis([-125.5, -122, 47, 49])

for sta in sta_df.index:
    x = sta_df.loc[sta, 'Longitude']
    y = sta_df.loc[sta, 'Latitude']
    
    ax.text(x, y, str(sta))

plt.show()