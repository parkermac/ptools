# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 14:13:03 2016

@author: PM5

This makes the file river_info.csv, and the individual river lon, lat: tracks.
"""

import os
import sys

alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()
dir0 = Ldir['parent']

import scipy.io as sio
import numpy as np
import pandas as pd

#%% make places for output
ri_name = 'pnw_all_2016_07'
ri_dir0 = dir0 + 'ptools_output/river/'
ri_dir = ri_dir0 + ri_name +'/'
rit_dir = ri_dir + 'tracks/'
Lfun.make_dir(ri_dir0, clean=False)
Lfun.make_dir(ri_dir, clean=True)
Lfun.make_dir(rit_dir, clean=True)

#%% load information from SNG
fn = (Ldir['data'] +
    'rivers/Info_curated/' +
    'SNG_riverInformation.mat')
aa = sio.loadmat(fn)
all_riv_lon = aa['lon']
all_riv_lat = aa['lat']
rr = aa['rivers']
NR, NC = rr.shape
name = []
depth = np.zeros(NC)
width = np.zeros(NC)
max_dist = np.zeros(NC)
ratio_sng = np.zeros(NC)
for ii in range(NC):
    a = rr[0, ii]
    # some renaming
    rn = a[0][0].lower()
    if rn == 'duwamish':
        rn = 'green'
    elif rn == 'hammahamma':
        rn = 'hamma'
    name.append(rn) # string
    #
    # save the river track information to a separate file
    lon = a[1].flatten() # ndarray
    lat = a[2].flatten() # ndarray
    df_tr = pd.DataFrame()
    df_tr['lon'] = lon
    df_tr['lat'] = lat
    df_tr.index.name = 'ind'
    fn_tr = rit_dir + rn + '.csv'
    df_tr.to_csv(fn_tr)
    #
    # save other information for use in a DataFrame
    depth[ii] = a[3].flatten()[0]
    width[ii] = a[4].flatten()[0]
    max_dist[ii] = a[8].flatten()[0]
    ratio_sng[ii] = a[9].flatten()[0]
    ii += 1

#%% initialize a DataFrame to organize the info
df = pd.DataFrame(index=name)

df['usgs'] = np.nan # final gage to use
df['usgs_sng'] = np.nan
df['usgs_ecy'] = np.nan # scale gage
df['usgs_nws'] = np.nan

df['ratio'] = np.nan
df['ratio_sng'] = ratio_sng
df['ratio_ecy'] = np.nan

df['nws'] = np.nan
df['ec'] = np.nan

df['depth'] = depth
df['width'] = width.astype('int')
df['max_dist'] = max_dist

#%% info from lists from Sarah

ec_code = ['08GB013', '08GA022', '08MF005', '08HB011', '08HD011', '08HB002',
           '08HA011', '08HB034', '08HA010', '08HB014', '08HC001']
ec_name = ['Clowhom', 'Squamish', 'Fraser', 'Tsolum', 'Oyster', 'Englishman',
           'Cowichan', 'Nanaimo', 'SanJuan', 'Sarita', 'Gold']
ec_dict = dict(zip(ec_name, ec_code))

for RN in ec_dict.keys():
    rn = RN.lower()
    df.loc[rn,'ec'] = ec_dict[RN]

usgs_code = ['12213100','12201500','12043300','12045500','12048000',
             '12054000','12061500','12079000','12089500','12101500',
             '12113000','12119000','12150800' ,'12167000','12200500',
             '12043000','12041200','12040500','12039500','12039005',
             '12031000','12013500','12010000','14246900','14301000',
             '14301500','14303600','14305500','14306500','14307620',
             '14321000','14325000']
usgs_name = ['Nooksack','Samish','Hoko','Elwha','Dungeness','Duckabush',
             'Skokomish','Deschutes','Nisqually','Puyallup','Green','Cedar',
             'Snohomish','Stillaguamish','Skagit','Calawah','Hoh','Queets',
             'Quinault','Humptulips','Chehalis','Willapa','Naselle',
             'Columbia','Nehalem','Wilson','Nestucca','Siletz','Alsea',
             'Siuslaw','Umpqua','Coquille']
             

usgs_dict = dict(zip(usgs_name, usgs_code))

# Add additional rivers by hand
usgs_dict['nf_skokomish'] = '12059500'
usgs_dict['sf_skokomish'] = '12060500'

for RN in usgs_dict.keys():
    rn = RN.lower()
    df.loc[rn,'usgs_sng'] = usgs_dict[RN]

#%% load other info

def get_nws_info(Ldir, riv_name):
    out_dict = dict()
    # This listing has all the stations with NWS Forecasts.
    fn = Ldir['data'] + 'rivers/Info_curated/USGS_NWS_Codes.csv'
    df = pd.read_csv(fn, index_col='Name')
    if riv_name in df.index:
        out_dict['usgs_nws'] = str(int(df.loc[riv_name, 'Station Number']))
        has_nws = df.loc[riv_name, 'NWS Forecast']
        if has_nws == 'YES':
            out_dict['nws'] = df.loc[riv_name, 'NWS ID']
    else:
        pass
    return out_dict

def get_ecy_info(Ldir, riv_name):
    out_dict = dict()
    # This listing, from Mohamedali et al. (2011) has a long
    # list of rivers coming into the Salish Sea, and then associates
    # each with a USGS gage and a scaling factor.
    # We cordinate this list using riv_name, and assume it is the same
    # as the lower case version of the first word in the index (Watershed Name).
    fn = Ldir['data'] + 'rivers/Info_curated/Ecology_Scale_Factors.csv'
    df = pd.read_csv(fn, index_col='Watershed Name')
    for wn in df.index:
        if wn.split()[0].lower() == riv_name:
            sg = df.loc[wn]['Scale Gage'].split()[-1]
            sg = sg[sg.find('(')+1 : sg.find(')')]
            out_dict['usgs_ecy'] = sg
            out_dict['ratio_ecy'] = df.loc[wn]['Scale Factor']
    return out_dict

# here we go through all the river names that are in the data frame
# and add other entries for them
for rn in df.index:
    try:
        out_dict_nws = get_nws_info(Ldir, rn)
        for item in out_dict_nws.keys():
            df.loc[rn, item] = out_dict_nws[item]
    except:
        pass
    try:
        out_dict_ecy = get_ecy_info(Ldir, rn)
        for item in out_dict_ecy.keys():
            if rn in ['fraser','green'] and item == 'usgs_ecy':
                pass
            else:
                df.loc[rn, item] = out_dict_ecy[item]
    except:
        pass

#%% decide which gages and ratios to use

# initialize with Sarah's list
df['usgs'] = df['usgs_sng']
df['ratio'] = df['ratio_sng']

# then replace any instances which have a scale gage from Ecology
for rn in df.index:
    if pd.notnull(df.loc[rn,'usgs_ecy']):
        df.loc[rn,'usgs'] = df.loc[rn,'usgs_ecy']
        df.loc[rn,'ratio'] = df.loc[rn,'ratio_ecy']
        print(rn + ': replacing usgs with usgs_ecy')

# and check if there is any mismatch between the final result
# and the usgs number from nws
for rn in df.index:
    if ( (df.loc[rn,'usgs'] != df.loc[rn,'usgs_nws'])
        and pd.notnull(df.loc[rn,'usgs_nws'])
        and pd.notnull(df.loc[rn,'usgs']) ):
        print(rn + ': ' + str(df.loc[rn,'usgs']) +
              '.ne.' + str(df.loc[rn,'usgs_nws']))

# clean up:

# drop rivers without any gauging station
rn_to_drop = []
for rn in df.index:
    if ( pd.isnull(df.loc[rn,'usgs']) and pd.isnull(df.loc[rn,'ec']) ):
        rn_to_drop.append(rn)
        print(' -- Dropping ' + rn)
df = df.drop(rn_to_drop)

# set any missing data to default values
# (written for the case of rivers I added by hand)
for rn in df.index:
    if pd.isnull(df.loc[rn,'ratio']):
        df.loc[rn,'ratio'] = 1.0
    if pd.isnull(df.loc[rn,'depth']):
        df.loc[rn,'depth'] = 5.0
    if pd.isnull(df.loc[rn,'width']):
        df.loc[rn,'width'] = 2.0
    if pd.isnull(df.loc[rn,'max_dist']):
        df.loc[rn,'max_dist'] = 0.1

#%% save to a csv file
fn_ri = ri_dir + 'river_info.csv'

df_final = df[['usgs', 'ec', 'nws', 'ratio', 'depth', 'width', 'max_dist']]
df_final.index.name = 'rname'

df_final.to_csv(fn_ri)

# here is how you would read it back into a DataFrame
# df1 = pd.read_csv(fn_ri, index_col='rname')


