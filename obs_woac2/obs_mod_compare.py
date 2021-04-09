"""
This reads in and processes the WOAC csv files from Simone Alin,
comparing them to LiveOcean extractions.

"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import netCDF4 as nc
from datetime import datetime
import matplotlib.pyplot as plt
import seawater as sw
from PyCO2SYS import CO2SYS

pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3', ex_name='lo8b')

# where are the model extractions
folder = 'woac'
mod_pth = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'cast' / folder

# where are the observations
obs_pth = Path(__file__).absolute().parent.parent.parent / 'ptools_data' / 'woac' / 'raw_2019_10'

# get the list of cruises and casts
sta_df_pth = Path(__file__).absolute().parent.parent.parent / 'ptools_output' / 'woac2'
sta_df = pd.read_pickle(sta_df_pth / 'sta_df.p')

plt.close('all')
cruises = sta_df['Cruise'].unique()
iii = 0
for cruise in cruises:
    print(cruise)
    
    fn = obs_pth / (cruise + '_data.csv')
    obs_df = pd.read_csv(fn, parse_dates=[['DATE_UTC','TIME_UTC']])
    # the index at this point is just row number
    
    # set known bad data to np.nan
    obs_df[obs_df=='nan nan'] = np.nan
    obs_df[obs_df==-999] = np.nan
    
    # drop rows or columns that have no good data at all, or that are empty
    obs_df = obs_df.dropna(how='all')
    
    # keep only selected columns
    obs_df = obs_df[['STATION_NO', 'CTDPRS_DBAR', 'CTDTMP_DEG_C_ITS90',
       'CTDSAL_PSS78', 'CTDOXY_UMOL_KG_ADJ', 'TA_UMOL_KG', 'DIC_UMOL_KG',
       'NITRATE_UMOL_KG', 'NITRITE_UMOL_KG', 'AMMONIA_UMOL_KG', 'AMMONIUM_UMOL_L',
       'PHOSPHATE_UMOL_KG', 'SILICATE_UMOL_KG']]
       
    # rename
        
    c_df = sta_df[sta_df['Cruise']==cruise]
    
    for sn in c_df.index:
        print(sn)
        
        # observations
        o_df = obs_df[obs_df['STATION_NO']==sn]
        lat = c_df.loc[sn, 'Latitude']
        O = dict()
        O['p'] = o_df['CTDPRS_DBAR'].values
        O['z'] = -sw.dpth(O['p'], lat)
        O['s'] = o_df['CTDSAL_PSS78'].values
        O['t'] = o_df['CTDTMP_DEG_C_ITS90'].values # this is in-situ
        O['th'] = sw.ptmp(O['s'], O['t'], O['p'])
        O['do'] = o_df['CTDOXY_UMOL_KG_ADJ'].values
        O['din'] = o_df[['NITRATE_UMOL_KG', 'NITRITE_UMOL_KG', 'AMMONIA_UMOL_KG',
                'AMMONIUM_UMOL_L']].sum(axis=1, min_count=1).values
        O['ta'] = o_df['TA_UMOL_KG'].values
        O['dic'] = o_df['DIC_UMOL_KG'].values
        
        # model
        mod_fn = mod_pth / (cruise + '_' + str(int(sn)) + '.nc')
        m_ds = nc.Dataset(mod_fn)
        #  variable names are:
        # AKs,salt,temp,NO3,phytoplankton,zooplankton,
        # detritus,Ldetritus,oxygen,alkalinity,TIC,h
        M = dict()
        M['s'] = m_ds['salt'][:].squeeze().data
        M['th'] = m_ds['temp'][:].squeeze().data
        M['z'] = m_ds['z_rho'][:].squeeze().data
        pdens = sw.dens(M['s'], M['th'], 0)
        M['do'] = m_ds['oxygen'][:].squeeze().data*1000/pdens
        M['din'] = m_ds['NO3'][:].squeeze().data*1000/pdens
        M['ta'] = m_ds['alkalinity'][:].squeeze().data*1000/pdens
        M['dic'] = m_ds['TIC'][:].squeeze().data*1000/pdens
        
        # pack into DataFrames
        OO = pd.DataFrame(data=O)
        OO = OO.sort_values('z')
        OO = OO.set_index(OO['z'])
        
        MMr = pd.DataFrame(data=M)
        MMr = MMr.sort_values('z')
        MMr = MMr.set_index(MMr['z'])
        
        # interpolate modeled data to observed depths
        MM = MMr.reindex(OO.index, method='nearest')
        
        # keep only selected variables
        v_list = ['z', 's', 'th', 'do', 'din', 'ta', 'dic']
        vv_list = ['s', 'th', 'do', 'din', 'ta', 'dic']
        OO = OO[v_list]
        MM = MM[v_list]
        MMr = MMr[v_list]
        
        # combine into a single DataFrame
        OM = OO.copy()
        for vn in vv_list:
            OM[vn+'_m'] = MM[vn]
        OM['Station'] = sn
        OM['Cruise'] = cruise
        
        if iii == 0:
            OMa = OM.copy()
        else:
            OMa = pd.concat((OMa, OM), ignore_index=True)
        
        iii += 1
        
        if False:
            fig = plt.figure(figsize=(20,6))
            vv_list = ['s', 'th', 'do', 'din', 'ta', 'dic']
            NV = len(vv_list)
            ii = 1
            for vn in vv_list:
                ax = fig.add_subplot(1, NV, ii)
                OO.plot(x=vn, y='z', title=vn, style='-ok', ax=ax, label='obs')
                MMr.plot(x=vn, y='z', style='oc', ax=ax, label='mod')
                MM.plot(x=vn, y='z', style='*r', ax=ax, label='modi')
                ax.set_ylim(-250,0)
                ax.grid(True)
                ii += 1

# calculate ph and arag (pH and aragonite saturation state)
CO2dict = CO2SYS(OMa['ta'].values, OMa['dic'].values, 1, 2,
    OMa['s'].values, OMa['th'].values, OMa['th'].values,
    0, 0, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)
CO2dict_m = CO2SYS(OMa['ta_m'].values, OMa['dic_m'].values, 1, 2,
    OMa['s_m'].values, OMa['th_m'].values, OMa['th_m'].values,
    0, 0, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)
OMa['ph'] = CO2dict['pHout']
OMa['ph_m'] = CO2dict_m['pHout']
OMa['arag'] = CO2dict['OmegaARout']
OMa['arag_m'] = CO2dict_m['OmegaARout']

OMa.to_pickle(sta_df_pth / 'OMa.p')
        
plt.show()
        

