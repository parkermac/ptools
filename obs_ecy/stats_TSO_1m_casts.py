"""
Calculates validation statistics and other quantities for temperature,
salinity and DO from Ecology and EC CTD casts.  Uses all depths in 1 m
increments.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart()
import zfun

# ***** User Edits

Ldir['gtagex'] = 'cas6_v3_lo8b'

ps_only = False

# ***** End User Edits

TSO_all = pd.DataFrame()
TSO_dict = dict()

for year in [2017, 2018, 2019]:
    
    # input location
    in_dir = Ldir['parent'] + 'ptools_output/ecology/'
    in_fn = 'TSO_1m_casts_' + Ldir['gtagex'] + '_' + str(year) + '.p'
    
    TSO = pd.read_pickle(in_dir + in_fn)
    
    # drop stations that are not in Puget Sound
    drop_list = ['GYS008', 'GYS016','WPA003', 'WPA004', 'WPA006', 'WPA007',
        'WPA008', 'WPA113', 'SOG27', 'SOG42']
    
    if ps_only:
        for sn in drop_list:
            TSO = TSO[TSO['Station'] != sn]
            
    TSO_all = pd.concat((TSO_all, TSO), sort=False)
    
    # statistics by year
    DS = TSO['Mod Salinity'] - TSO['Salinity']
    DT = TSO['Mod Temp. (deg C)'] - TSO['Temp. (deg C)']
    DO = TSO['Mod DO (mg L-1)'] - TSO['DO (mg L-1)']

    print('\n*** ' + str(year) + ' ***')
    if ps_only:
        print('&&& Using only Puget Sound Stations &&&')
    print('Salinity Bias = %0.2f, RMSE = %0.2f' % (DS.mean(), np.sqrt((DS**2).mean())) )
    print('Mean Salinity = %0.2f, mean model salinity = %0.2f' % 
        (TSO['Salinity'].mean(), TSO['Mod Salinity'].mean()))
    print('Temp. (deg C) Bias = %0.2f, RMSE = %0.2f' % (DT.mean(), np.sqrt((DT**2).mean())) )
    print('DO (mg L-1) Bias = %0.2f, RMSE = %0.2f' % (DO.mean(), np.sqrt((DO**2).mean())) )

# statistics all years
DS = TSO_all['Mod Salinity'] - TSO_all['Salinity']
DT = TSO_all['Mod Temp. (deg C)'] - TSO_all['Temp. (deg C)']
DO = TSO_all['Mod DO (mg L-1)'] - TSO_all['DO (mg L-1)']

print('\n*** ALL YEARS ***')
if ps_only:
    print('&&& Using only Puget Sound Stations &&&')
print('Salinity Bias = %0.2f, RMSE = %0.2f' % (DS.mean(), np.sqrt((DS**2).mean())) )
print('Temp. (deg C) Bias = %0.2f, RMSE = %0.2f' % (DT.mean(), np.sqrt((DT**2).mean())) )
print('DO (mg L-1) Bias = %0.2f, RMSE = %0.2f' % (DO.mean(), np.sqrt((DO**2).mean())) )

