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
from datetime import datetime, timedelta
import seawater

# where the data is, and where results will be stored
dir0 = '../../ptools_data/woac/'

a_raw = pd.read_excel(dir0 + 'raw/WOAC_data_9-7-2018_dataToParker.xlsx',
    parse_dates = ['Date_collected', 'Time_collected'])
    
# keep only selected columns
a = a_raw[['CRUISE_ID', 'Date_collected', 'Time_collected',
       'LONGITUDE_DEC', 'LATITUDE_DEC', 'STATION_NO',
       'CTDPRS_DBAR', 'CTDTMP_DEG_C_ITS90', 'CTDSAL_PSS78',
       'SIGMATHETA_KG_M3',
       'CTDOXY_UMOL_KG_ADJ', 'OXYGEN_FLAG_W',
       'NITRATE_UMOL_L', 'NITRITE_UMOL_L', 'AMMONIUM_UMOL_L',
       'CHLA avg (ug/l)',
       'TA_1_umol_kg', 'DIC_1_umol_kg',
       'TA_1_QC ', 'DIC_1_QC',
       'TA_2_umol_kg', 'DIC_2_umol_kg',
       'TA_2_QC ', 'DIC_2_QC',
       'Ph Total in situ', 'pCO2 µatm', 'CO2 µmol/kg',
       'HCO3- µmol/kg', 'CO3-- µmol/kg',
       'Omega Ca', 'Omega Ar']]
       
# rename some columns
a = a.rename(columns={'Date_collected':'Date', 'Time_collected':'Time',
    'STATION_NO':'Station', 'CRUISE_ID':'Cruise',
    'LONGITUDE_DEC':'Longitude', 'LATITUDE_DEC':'Latitude',
    'CTDPRS_DBAR':'Pressure (dbar)', 'CTDTMP_DEG_C_ITS90':'Temp. (deg C)',
    'CTDSAL_PSS78':'Salinity','SIGMATHETA_KG_M3':'Sigma (kg m-3)',
    'CTDOXY_UMOL_KG_ADJ':'DO (umol kg-1)', 'OXYGEN_FLAG_W':'DO Flag',
    'NITRATE_UMOL_L':'NO3 (uM)', 'NITRITE_UMOL_L':'NO2 (uM)',
    'AMMONIUM_UMOL_L':'NH4 (uM)',
    'CHLA avg (ug/l)':'Chl (ug L-1)',
    'TA_1_umol_kg':'TA1 (umol kg-1)', 'DIC_1_umol_kg':'DIC1 (umol kg-1)',
    'TA_1_QC ':'TA1 Flag', 'DIC_1_QC':'DIC1 Flag',
    'TA_2_umol_kg':'TA2 (umol kg-1)', 'DIC_2_umol_kg':'DIC2 (umol kg-1)',
    'TA_2_QC ':'TA2 Flag', 'DIC_2_QC':'DIC2 Flag',
    'Ph Total in situ':'pH', 'pCO2 µatm':'pCO2 (uatm)',
    'CO2 µmol/kg':'CO2 (umol kg-1)', 'HCO3- µmol/kg':'HCO3- (umol kg-1)',
    'CO3-- µmol/kg':'CO3-- (umol kg-1)'
    })
    
# nans for bad data
a[a==-999.] = np.nan

# convert units:
a['DO (uM)'] = a['DO (umol kg-1)'] * (1000 + a['Sigma (kg m-3)'])/1000
a.pop('DO (umol kg-1)')

# derived quantities
a['Z (m)'] = -seawater.dpth(a['Pressure (dbar)'].values, a['Latitude'].values)

a['Datetime'] = a['Date']  + pd.to_timedelta(a['Time'])
a.pop('Date')
a.pop('Time')

# number the casts
castnum = 0
sta0 = a.loc[0,'Station']
for ii in a.index:
    sta = a.loc[ii,'Station']
    if sta == sta0:
        pass
    else:
        castnum += 1
        sta0 = a.loc[ii,'Station']
    a.loc[ii,'castnum'] = castnum

# save result
a.to_pickle(dir0 + 'Casts_2017.p')

# then associate casts with station number and mean lat, lon, time
sta_df = pd.DataFrame(columns = ['Station','Cruise','Latitude', 'Longitude','Datetime'])
sta_df.index.name = 'castnum'
for cast in set(a['castnum']):
    aa = a[a['castnum'] == cast]
    nn = len(aa)
    sta_df.loc[cast,'Longitude'] = aa.loc[:,'Longitude'].mean()
    sta_df.loc[cast,'Latitude'] = aa.loc[:,'Latitude'].mean()
    sta_df.loc[cast,'Station'] = aa.loc[aa.index[0],'Station']
    sta_df.loc[cast,'Cruise'] = aa.loc[aa.index[0],'Cruise']
    sta_df.loc[cast,'Datetime'] = aa.loc[aa.index[0] + int(nn/2),'Datetime']
# save sta_df for later use
sta_df.to_pickle(dir0 + 'sta_df.p')

