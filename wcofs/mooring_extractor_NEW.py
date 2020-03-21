"""
Extract multiple mooring-like records from WCOFS.

The output is structured to be like LiveOcean/x_moor/layer_extractor.py,
and so can be plotted with x_moor/plot_mooring.py.

This is designed to work with the files downloaded areound March 2020.

Because these have a lat-lon grid provided we will just use nearest neighbor
interpolation to find nearby points, and not worry about the distinction
between the various C-grids.

NOTE: in three instances the file size is larger, and there are two time levels
instead of one.  This corresponds to times when the day before is missing, so
I assume that these contain the "missing".

Oddly, the timestamp for nos.wcofs.avg.nowcast.20190122.t03z.nc implies:
datetime(2016,1,1) + timedelta(days=96476400/86400)
Out[13]: datetime.datetime(2019, 1, 21, 15, 0)
that the time is 3 PM from the day before...?  Perhaps this is a side effect of the
time averaging, like it is the start of the averaging time.
"""

testing = True

# setup
import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
import numpy as np
from datetime import datetime, timedelta
start_time = datetime.now()
import zfun
import zrfun

from scipy.spatial import cKDTree
import netCDF4 as nc

import wcofs_fun as wfun
from importlib import reload
reload(wfun)

# load LO mooring functions and lists
sys.path.append(os.path.abspath('../../LiveOcean/x_moor'))
import moor_fun as mfun
reload(mfun)
import moor_lists as ml
reload(ml)

# get the station locations
job_name = 'comt3_2014_offshore' # used in ml.get_sta_dict()
sta_dict, v2_list, v3_list_rho, v3_list_w = ml.get_sta_dict(job_name)

# form mooring location vectors
xs_list = []
ys_list = []
for sn in sta_dict.keys():
    xs_list.append(sta_dict[sn][0])
    ys_list.append(sta_dict[sn][1])
xs = np.array(xs_list)
ys = np.array(ys_list)
xys = np.array((xs,ys)).T
NS = len(sta_dict)

# get list of files to extract from
Ldir = Lfun.Lstart()
Ldir['gtagex'] = 'wcofs_avg_now'
in_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/'
fn_list_raw = os.listdir(in_dir)
fn_list = [(in_dir + ff) for ff in fn_list_raw if 'nos.wcofs.avg.nowcast' in ff]
fn_list.sort()
NT = len(fn_list)

# get grid info
G = wfun.get_grid_info(fn_list[0])

# create the nearest neighbor Tree object
# (the speed of the method relies on doing this)
mr = G['mask_rho'] == 1. # water points
mu = G['mask_u'] == 1. # water points
mv = G['mask_v'] == 1. # water points

xxr = G['lon_rho'][mr]; yyr = G['lat_rho'][mr]
xxu = G['lon_u'][mu]; yyu = G['lat_u'][mu]
xxv = G['lon_v'][mv]; yyv = G['lat_v'][mv]

XYr = np.array((xxr,yyr)).T
XYTr = cKDTree(XYr)

XYu = np.array((xxu,yyu)).T
XYTu = cKDTree(XYu)

XYv = np.array((xxv,yyv)).T
XYTv = cKDTree(XYv)

# get nearest lon, lat on each grid, as a test of the method
xrs = xxr[XYTr.query(xys, n_jobs=-1)[1]]
yrs = yyr[XYTr.query(xys, n_jobs=-1)[1]]
xus = xxu[XYTu.query(xys, n_jobs=-1)[1]]
yus = yyu[XYTu.query(xys, n_jobs=-1)[1]]
xvs = xxv[XYTv.query(xys, n_jobs=-1)[1]]
yvs = yyv[XYTv.query(xys, n_jobs=-1)[1]]
if testing == True:
    sys.path.append(os.path.abspath('../../LiveOcean/plotting'))
    import pfun
    import matplotlib.pyplot as plt
    plt.close('all')
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    # plot station locations
    ax.plot(xs, ys, 'ok')
    # plot nearest unmasked points on all three grids
    ax.plot(xrs, yrs, '+r')
    ax.plot(xus, yus, '+b')
    ax.plot(xvs, yvs, '+g')
    # and plot the underlying grid
    alpha = .2
    ax.plot(xxr, yyr, 'or', alpha=alpha)
    ax.plot(xxu, yyu, 'ob', alpha=alpha)
    ax.plot(xxv, yyv, 'og', alpha=alpha)
    # formatting
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-126, -123, 44, 49])
    plt.show()

fn = fn_list[0]
ds = nc.Dataset(fn)

v = ds['salt'][0,:,:,:]
N = v.shape[0]

VS = np.nan * np.ones((N,NS))
for n in range(N):
    vv = v[n,:,:][mr].data
    VS[n,:] = vv[XYTr.query(xys, n_jobs=-1)[1]]


