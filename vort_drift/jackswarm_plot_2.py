"""
More Jackswarm plotting.
"""
# setup
import os; import sys
alp = os.path.abspath('/Users/pm7/Documents/LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
from importlib import reload
import zfun; reload(zfun)
import matfun; reload(matfun)

import numpy as np
import matplotlib.pyplot as plt

def jack_to_df(jnum):
    # Function to load Jack data into a data frame
    dir0 = '/Users/pm7/Documents/tools_data/vortex_drifter/2014_11_JackSwarm_Data/'
    fn = 'JACK' + str(jnum) + '.CSV'
    # load a jack data frame, and make a datetime column
    import pandas as pd
    jdf = pd.read_csv(dir0 + fn, parse_dates = [['Date',' Time']])
    # remove spaces from column headings (how does a lambda function work?)
    cols = jdf.columns
    cols = cols.map(lambda x: x.replace(' ', '') if isinstance(x, (str, unicode)) else x)
    jdf.columns = cols    
    # and specify the index column
    jdf1 = jdf.set_index('Date_Time')        
    return jdf1

# list of jacks and their colors
jlist = [108, 205, 317, 405, 511]
clist = ['r','b','g','m','c']
cdict = dict(zip(jlist,clist))

# put all the data frames in a dict
Jdf = dict()
for jnum in jlist:
    Jdf[jnum] = jack_to_df(jnum)

# plotting
plt.close()
fig = plt.figure(figsize=(15,8))
ax = fig.add_subplot(111)

for jnum in jlist:
    jdf = Jdf[jnum]
    jdf1 = jdf[jdf['Mode']==2] # pull out times in drift mode only
    dt = jdf1.index.values
    # get a single starting time, from the first jack in the list
    if jnum == jlist[0]:
        dt0 = dt[0]
    head = jdf1['Heading'].values
    # center point for the x-y grid
    lon0 = np.deg2rad(-122.9)
    lat0 = np.deg2rad(47.16)
    RE = 6371e3 # Earth radius (m)
    x = RE*(np.deg2rad(jdf1['Longitude'].values) - lon0)*np.cos(lat0)
    y = RE*(np.deg2rad(jdf1['Latitude'].values) - lat0)
    # make a time axis out of the datetime format index
    ttd = (dt-dt0) # array of timedeltas
    tf = ttd.astype('float') # convert to array of floats
    ts = tf/1e9 # convert to seconds (timedeltas are in nanoseconds (1e-9 sec))
    # vorticity
    dts = np.diff(ts)
    dth = np.diff(head)        
    # fix places in the the dth vector that pass over the 0-360 jump
    dth0 = 180; # min size of angular jump to consider
    maskp = dth>dth0; maskn = dth<-dth0        
    dthp = dth; dthp[maskp] = dth[maskp] - 360
    dthpn = dthp; dthpn[maskn] = 360 + dthp[maskn]
    vort = -2*(np.pi/180) * dthpn / dts;        
    vort[dts > 5] = np.nan # mask gaps between drifts

    ax.plot(x, y, '.', color=cdict[jnum])
    
plt.show()