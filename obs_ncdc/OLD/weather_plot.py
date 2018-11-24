"""
This plots met data.  It is aimed at exploring decade-to-century
variation in conditions in Puget Sound.
"""
import os
import shutil
import pandas as pd

data_dir0 = '../../ptools_data/ncdc/'
fn = data_dir0 + 'seatac_1948-2015.csv'

def make_dir(dirname, clean=False):
    if clean == True:
        shutil.rmtree(dirname, ignore_errors=True)
        os.mkdir(dirname)
    else:
        try:
            os.mkdir(dirname)
        except OSError:
            pass

do_save = True
if do_save:
    out_dir00 = '../../ptools_output/'
    out_dir0 = out_dir00 + 'ncdc/'
    make_dir(out_dir00)
    make_dir(out_dir0)
    out_fn = out_dir0 + 'seatac_monthly_1948-2015.csv'
    out_fnp = out_dir0 + 'seatac_monthly_1948-2015.p'

# read csv data into a data frame
df = pd.read_csv(fn, index_col = 'DATE', parse_dates = ['DATE'])

import numpy as np
df[df == -9999] = np.nan

cols = list(df.columns)
sta_name = df.loc[0,'STATION_NAME']
sta_id = df.loc[0,'STATION']

# make a list of tuples with:
# (variable name, scaling factor, label text, low, high)

vtup_list = [
    ('TMAX',0.1,'Tmax (degC)',0,30),
    ('TMIN',0.1,'Tmin (degC)',-10,20),
    ('PRCP',0.1,'Precip (mm/day)',0,10),
    ('ACMH',1.0,'Clouds (obs) %',0,100),
    ('ACSH',1.0,'Daytime Clouds (obs) %',0,100),
    ('PSUN',1.0,'% of Possible Sunshine',0,100),
    ('TSUN',1/60.,'Daily Sunshine (hours)',0,12),
    ('AWND',0.1,'Average Daily Windspeed (m/s)',0,5),
    ('WDF2',1.0,'Dir of 2-minute Wind (deg)',0,360)
    ]

vname_list = []
for tup in vtup_list:
    vname_list.append(tup[0])

from datetime import datetime

df = df[vname_list]

for tup in vtup_list:
    vn = tup[0]
    df[vn] = df[vn]*tup[1]

# get monthly and annual means
dfm = df.resample('M', label='left', loffset='15d').mean()
dfm.TSUN[dfm.TSUN==0] = np.nan
dfy = df.resample('A', label='left', loffset='6m').mean()
dfyj = df.resample('A-JUN', label='left', loffset='6m').mean()

# # climatology
# dfma = dfm.copy()
# for mm in range(1,13):
#     im = dfm.index.month == mm
#     dfma[im] = dfm[im] - dfm[im].mean()

import matplotlib.pyplot as plt
plt.close()

NR = 3; NC = 3
fig = plt.figure(figsize=(20,12))

xlimtup = (datetime(1945,1,1),datetime(2017,1,1))

cc = 1
for tup in vtup_list:
    vn = tup[0]
    lab = tup[2]
    axx = fig.add_subplot(NR,NC,cc)
    axx.set_ylabel(lab)
    axx.set_ylim(tup[3],tup[4])
    dfm.plot(ax = axx, y = vn,
        legend=False, color='b',xlim=xlimtup)
    if vn == 'PRCP':
        dfyj.plot(ax = axx, y = vn, legend=False, color='r',linewidth=3,xlim=xlimtup)
    else:
        dfy.plot(ax = axx, y = vn, legend=False, color='r',linewidth=3,xlim=xlimtup)
    cc += 1

plt.show()

if do_save:
    dfm.to_csv(out_fn)
    dfm.to_pickle(out_fnp)
