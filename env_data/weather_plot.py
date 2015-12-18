"""
This plots met data.  It is aimed at exploring decade-to-century
variation in conditions in Puget Sound.
"""

dir0 = '/Users/PM5/Documents/tools_data/obs_data/met/'
fn = dir0 + 'seatac_1948-2015_v2.csv'
# read csv data into a data frame
import pandas as pd    
df = pd.read_csv(fn, index_col = 'DATE', parse_dates = ['DATE'])

import numpy as np
df[df == -9999] = np.nan

cols = list(df.columns)
sta_name = df.ix[0,'STATION_NAME']
sta_id = df.ix[0,'STATION']

# make a list of tuples with:
# (variable name, scaling factor, label text, low, high)


vtup_list = [
    ('TMAX',0.1,'Tmax (degC)',0,30),
    ('TMIN',0.1,'Tmin (degC)',-10,20),
    ('PRCP',0.1,'Precip (mm/day)',0,10),
    ('ACMH',1.0,'Clouds (obs) %',0,100),
    #('ACSH',1.0,'Daytime Clouds (obs) %',0,100),
    ('PSUN',1.0,'% of Possible Sunshine',0,100),
    #('TSUN',1/60.,'Daily Sunshine (hours)',0,12),
    #('AWND',0.1,'Average Daily Windspeed (m/s)',0,5),
    #('WDF2',1.0,'Dir of 2-minute Wind (deg)',0,360)
    ]

vname_list = []
for tup in vtup_list:
    vname_list.append(tup[0])

from datetime import datetime
df = df[df.index<datetime(2015,1,1)]

df = df[vname_list]

for tup in vtup_list:
    vn = tup[0]
    df[vn] = df[vn]*tup[1]

# get monthly and annual means
dfm = df.resample('M', how='mean', label='left', loffset='15d')
dfy = df.resample('A', how='mean', label='left', loffset='6m')
dfyj = df.resample('A-JUN', how='mean', label='left', loffset='6m')

# climatology?
dfma = dfm.copy()
for mm in range(1,13):
    im = dfm.index.month == mm
    dfma[im] = dfm[im] - dfm[im].mean()

import matplotlib.pyplot as plt
plt.close()

NR = 3; NC = 2
fig = plt.figure(figsize=(25,15))

xlimtup = (datetime(1945,1,1),datetime(2015,1,1))

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
    #dfma.plot(ax = axx, y = vn, legend=False, color='g')
    cc += 1

plt.show()
