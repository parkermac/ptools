"""
Plots mooring records to compare with IRIS record.
"""

# setup
import os
import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
from importlib import reload
import Lfun
reload(Lfun)
import numpy as np
import zfun
reload(zfun) 
import matfun
reload(matfun)
from datetime import datetime, timedelta
import pickle

# Get the model mooring record
Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'moor/'
fn = 'cascadia1_base_lo1_FN03C_low_pass_2013.01.02_2015.11.01.p'
V, v1_list, v2_list, v3_list, G = pickle.load( open( indir + fn, 'rb' ) )
# convert model time to datetime
mdt = []
for mt in V['ocean_time']:
    mdt.append(datetime(1970,1,1,0,0) + timedelta(seconds=mt))
mT = V['temp'][0,:]

# Get IRIS data
indir_iris = '/Users/PM5/Documents/tools_data/obs_data/mooring/IRIS/'
iris_sta = 'FN03C'
fn_iris = 'FN03C.15min.T.mat'
iris = matfun.loadmat(indir_iris + fn_iris)
# station information
iris_meta = matfun.loadmat(indir_iris + '7D_metadata.mat')
isn = iris_meta['sta_name']
ii = np.where(isn == iris_sta)[0][0]
ilon = iris_meta['lon'][ii]
ilat = iris_meta['lat'][ii]
# convert matlab datenum to python datetime
idn = iris['t0'] + iris['t']
idt = []
for dn in idn:
    idt.append( datetime.fromordinal(dn.astype(int)) +
        timedelta(days = np.mod(dn,1)) - timedelta(days = 366) )
iT = iris['T']

# plottting       
import matplotlib.pyplot as plt

plt.close()
NR = 1; NC = 1
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(15,8), squeeze=False)
ax = axes[0,0]
line1, = ax.plot(idt, iT, '-r', label='hi')
line2, = ax.plot(mdt, mT, '-b', label='hi')
line1.set_label('IRIS')
line2.set_label('LiveOcean')
ax.legend()
ax.set_ylim(6,14)
ax.grid()
ax.set_title(iris_sta)
ax.set_xlabel('Date')
ax.set_ylabel('Bottom Temperature $^{\circ}C$')
plt.show()


