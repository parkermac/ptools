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

testing = False

# setup
import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart()
Ldir['gtagex'] = 'wcofs_avg_now'

import numpy as np
from datetime import datetime, timedelta
start_time = datetime.now()
import zfun
import zrfun

from scipy.spatial import cKDTree
import netCDF4 as nc
import pandas as pd

import wcofs_fun as wfun
from importlib import reload
reload(wfun)

# load LO mooring functions and lists
sys.path.append(os.path.abspath('../../LiveOcean/x_moor'))
import moor_fun as mfun
reload(mfun)
import moor_lists as ml
reload(ml)

# get list of files to extract from
in_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/'
fn_list_raw = os.listdir(in_dir)
fn_list = [(in_dir + ff) for ff in fn_list_raw if 'nos.wcofs.avg.nowcast' in ff]
fn_list.sort()
NT = len(fn_list)

# get the station locations
job_name = 'comt3_2014_offshore'
sta_dict, v2_list, v3_list_rho, v3_list_w = ml.get_sta_dict(job_name)
# convert sta_dict to a DataFrame
sta_df = pd.DataFrame(0, index=sta_dict.keys(), columns=['Lon', 'Lat'])
for sn in sta_dict.keys():
    sta_df.loc[sn, 'Lon'] = sta_dict[sn][0]
    sta_df.loc[sn, 'Lat'] = sta_dict[sn][1]

# name output files
out_fn_dict = dict()
for sn in sta_df.index:
    # make sure the output directory exists
    outdir0 = Ldir['LOo'] + 'moor/'
    Lfun.make_dir(outdir0)
    outdir = outdir0 + Ldir['gtagex'] + '/'
    Lfun.make_dir(outdir)
    # name output file
    out_fn = (outdir + sn + '.nc')
    # get rid of the old version, if it exists
    try:
        os.remove(out_fn)
    except OSError:
        pass
    out_fn_dict[sn] = out_fn


# form mooring location vectors
xs = sta_df['Lon'].to_numpy()
ys = sta_df['Lat'].to_numpy()
xys = np.array((xs,ys)).T
NS = len(sta_dict)

# get grid info
G, S = wfun.get_grid_info(fn_list[0])
N = S['N']

# create the nearest neighbor Tree objects
# masks
mr = G['mask_rho'] == 1. # water points
mu = G['mask_u'] == 1. # water points
mv = G['mask_v'] == 1. # water points
# water points in all three grids
xxr = G['lon_rho'][mr]; yyr = G['lat_rho'][mr]
xxu = G['lon_u'][mu]; yyu = G['lat_u'][mu]
xxv = G['lon_v'][mv]; yyv = G['lat_v'][mv]
# trees
XYr = np.array((xxr,yyr)).T
XYTr = cKDTree(XYr)
XYu = np.array((xxu,yyu)).T
XYTu = cKDTree(XYu)
XYv = np.array((xxv,yyv)).T
XYTv = cKDTree(XYv)

# get nearest lon, lat on each grid, as well as h and angle
for vn in ['h', 'lon_rho', 'lat_rho','angle']:
    sta_df.loc[:,vn] = G[vn][mr][XYTr.query(xys, n_jobs=-1)[1]]
for vn in ['lon_u', 'lat_u']:
    sta_df.loc[:,vn] = G[vn][mu][XYTu.query(xys, n_jobs=-1)[1]]
for vn in ['lon_v', 'lat_v']:
    sta_df.loc[:,vn] = G[vn][mv][XYTv.query(xys, n_jobs=-1)[1]]

if testing == True:
    sys.path.append(os.path.abspath('../../LiveOcean/plotting'))
    import pfun
    import matplotlib.pyplot as plt
    plt.close('all')
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    
    # plot station locations
    sta_df.plot(x='Lon', y='Lat', style='*k', ax=ax)
    sta_df.plot(x='lon_rho', y='lat_rho', style='or', ax=ax)
    sta_df.plot(x='lon_u', y='lat_u', style='>g', ax=ax)
    sta_df.plot(x='lon_v', y='lat_v', style='^b', ax=ax)
    for sn in sta_df.index:
        ax.text(sta_df.loc[sn,'Lon'], sta_df.loc[sn,'Lat'], sn)
    
    # formatting
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-126, -123, 44, 49])
        
    plt.show()

# EXAMPLE Getting a 3D variable at all stations, at one time
# fn = fn_list[0]
# ds = nc.Dataset(fn)
# v = ds['salt'][0,:,:,:]
# VS = np.nan * np.ones((N,NS))
# for n in range(N):
#     vv = v[n,:,:][mr].data
#     VS[n,:] = vv[XYTr.query(xys, n_jobs=-1)[1]]
    
# ================================

# creat lists of variables to get
v0_list = ['h', 'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v', 'angle']
v1_list = ['ocean_time']
v2_list = ['zeta']
v3_list_rho = ['salt', 'temp', 'u', 'v']
v3_list_w = []

# start by getting info from the first file
ds = nc.Dataset(fn_list[0])

# generate dictionaries of long names and units
V_long_name = dict()
V_units = dict()
v_all_list = v1_list + v2_list + v3_list_rho + v3_list_w
for vv in v0_list:
    try:
        V_long_name[vv] = ds.variables[vv].long_name
    except:
        V_long_name[vv] = ''
    try:
        V_units[vv] = ds.variables[vv].units
    except:
        V_units[vv] = ''
for vv in v_all_list:
    try:
        V_long_name[vv] = ds.variables[vv].long_name
    except:
        V_long_name[vv] = ''
    try:
        V_units[vv] = ds.variables[vv].units
    except:
        V_units[vv] = ''
ds.close()

# start netcdf files
for sn in sta_df.index:
    out_fn = out_fn_dict[sn]
    mfun.start_netcdf(out_fn, S['N'], NT, v0_list, v1_list, v2_list,
        v3_list_rho, v3_list_w, V_long_name, V_units)

# save static variables
for sn in sta_df.index:
    out_fn = out_fn_dict[sn]
    print('Initializing ' + out_fn)
    foo = nc.Dataset(out_fn, 'a')
    for vv in v0_list:
        foo[vv][:] =  sta_df.loc[sn,vv]
    foo.close()
    
# EXTRACT TIME-DEPENDENT FIELDS
count = 0
for fn in fn_list:
    ds = nc.Dataset(fn)
    print(' working on %d of %d' % (count, NT))
    sys.stdout.flush()
    
    v = ds['zeta'][0,:,:][mr].data
    zeta_s = v[XYTr.query(xys, n_jobs=-1)[1]]
        
    v = ds['salt'][0,:,:,:]
    VS = np.nan * np.ones((N,NS))
    for n in range(N):
        vv = v[n,:,:][mr].data
        VS[n,:] = vv[XYTr.query(xys, n_jobs=-1)[1]]
        
    v = ds['temp'][0,:,:,:]
    VT = np.nan * np.ones((N,NS))
    for n in range(N):
        vv = v[n,:,:][mr].data
        VT[n,:] = vv[XYTr.query(xys, n_jobs=-1)[1]]
    
    v = ds['u'][0,:,:,:]
    VU = np.nan * np.ones((N,NS))
    for n in range(N):
        vv = v[n,:,:][mu].data
        VU[n,:] = vv[XYTu.query(xys, n_jobs=-1)[1]]

    v = ds['v'][0,:,:,:]
    VV = np.nan * np.ones((N,NS))
    for n in range(N):
        vv = v[n,:,:][mv].data
        VV[n,:] = vv[XYTv.query(xys, n_jobs=-1)[1]]
    
    scount = 0
    for sn in sta_df.index:
        out_fn = out_fn_dict[sn]
        foo = nc.Dataset(out_fn, 'a')
        v = ds['ocean_time'][:].squeeze()
        foo['ocean_time'][count] = v
        foo['zeta'][count] = zeta_s[scount]
        foo['salt'][count,:] = VS[:,scount]
        foo['temp'][count,:] = VT[:,scount]
        foo['u'][count,:] = VU[:,scount]
        foo['v'][count,:] = VV[:,scount]
        foo.close()
    count += 1

    ds.close()
# END OF EXTRACTING TIME-DEPENDENT FIELDS

# create z_rho and z_w (has to be done after we have zeta)
# also rotate velocity to be u = E-W, v = N-S
for sn in sta_df.index:
    out_fn = out_fn_dict[sn]
    foo = nc.Dataset(out_fn, 'a')

    zeta = foo['zeta'][:].squeeze()
    hh = foo['h'][:] * np.ones_like(zeta)
    z_rho, z_w = zrfun.get_z(hh, zeta, S)
    
    v_var = foo.createVariable('z_rho', float, ('ocean_time','s_rho'))
    v_var.long_name = 'z on rho points (positive up)'
    v_var.units = 'm'
    v_var[:] = z_rho.T
    
    v_var = foo.createVariable('z_w', float, ('ocean_time','s_w'))
    v_var.long_name = 'z on w points (positive up)'
    v_var.units = 'm'
    v_var[:] = z_w.T
    
    u = foo['u'][:]
    v = foo['v'][:]
    angle = foo['angle'][:] # scalar
    # long_name: angle between xi axis and east
    # units: degrees
    # we use the minus sign to get back to Earth coordinates
    ca = np.cos(-np.pi*angle/180)
    sa = np.sin(-np.pi*angle/180)
    uu = ca*u + sa*v
    vv = ca*v - sa*u
    foo['u'][:] = uu
    foo['v'][:] = vv
    
    foo.close()



