"""
Extract multiple mooring-like records from WCOFS.

The output is structured to be like LiveOcean/x_moor/layer_extractor.py,
and so can be plotted with x_moor/plot_mooring.py.
"""

testing = False

# setup
import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import numpy as np
from datetime import datetime, timedelta
start_time = datetime.now()
import zfun
import zrfun
import netCDF4 as nc

pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

import wcofs_fun as wfun
from importlib import reload
reload(wfun)

# load LO mooring functions and lists
pth = os.path.abspath('../../LiveOcean/x_moor')
if pth not in sys.path:
    sys.path.append(pth)
import moor_fun as mfun
reload(mfun)
import moor_lists as ml
reload(ml)

# get the station locations
job_name = 'comt3_2014_offshore' # used in ml.get_sta_dict()
sta_dict, v2_list, v3_list_rho, v3_list_w = ml.get_sta_dict(job_name)

# get list of files to extract from
Ldir = Lfun.Lstart()
Ldir['gtagex'] = 'WCOFS_avg_Exp37'
in_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/'
fn_list_raw = os.listdir(in_dir)
fn_list = [(in_dir + ff) for ff in fn_list_raw if 'zuvts_ORWA_Parker_Exp37' in ff]
fn_list.sort()
NT = len(fn_list)

if testing:
    fn_list = [fn_list[0]]
    sta_list = ['CE015']
else:
    sta_list = list(sta_dict.keys())
    
# name output files
out_fn_dict = dict()
for sta_name in sta_list:
    # make sure the output directory exists
    outdir0 = Ldir['LOo'] + 'moor/'
    Lfun.make_dir(outdir0)
    outdir = outdir0 + Ldir['gtagex'] + '/'
    Lfun.make_dir(outdir)
    # name output file
    out_fn = (outdir + sta_name + '.nc')
    # get rid of the old version, if it exists
    try:
        os.remove(out_fn)
    except OSError:
        pass
    out_fn_dict[sta_name] = out_fn

# vertical grid information
S_info_dict = {'THETA_S':8,
        'THETA_B':3,
        'TCLINE':50,
        'N':40,
        'VTRANSFORM':2,
        'VSTRETCHING':4}
S = zrfun.get_S(S_info_dict)

# get grid info (the whole WCOFS domain)
grid_fn = in_dir + 'grd_wcofs_large_visc200.nc'
G = wfun.get_grid_info(grid_fn)
# specify and save grid info just for the cut-out that Alex made for me
cut_j0 = 1030-1
cut_j1 = 1710
cut_i0 = 375-1
cut_i1 = 696
h = G['h'][cut_j0:cut_j1, cut_i0:cut_i1]
angle = G['angle'][cut_j0:cut_j1, cut_i0:cut_i1]
lon_rho = G['lon_rho'][cut_j0:cut_j1, cut_i0:cut_i1]
lat_rho = G['lat_rho'][cut_j0:cut_j1, cut_i0:cut_i1]
mask_rho = G['mask_rho'][cut_j0:cut_j1, cut_i0:cut_i1]
lon_u = G['lon_u'][cut_j0:cut_j1, cut_i0:cut_i1-1]
lat_u = G['lat_u'][cut_j0:cut_j1, cut_i0:cut_i1-1]
mask_u = G['mask_u'][cut_j0:cut_j1, cut_i0:cut_i1-1]
lon_v = G['lon_v'][cut_j0:cut_j1-1, cut_i0:cut_i1]
lat_v = G['lat_v'][cut_j0:cut_j1-1, cut_i0:cut_i1]
mask_v = G['mask_v'][cut_j0:cut_j1-1, cut_i0:cut_i1]
lon_psi = G['lon_psi'][cut_j0:cut_j1-1, cut_i0:cut_i1-1]
lat_psi = G['lat_psi'][cut_j0:cut_j1-1, cut_i0:cut_i1-1]
Gcut = dict()
Gcut['h'] = h
Gcut['angle'] = angle
Gcut['lon_rho'] = lon_rho
Gcut['lat_rho'] = lat_rho
Gcut['mask_rho'] = mask_rho
Gcut['lon_u'] = lon_u
Gcut['lat_u'] = lat_u
Gcut['mask_u'] = mask_u
Gcut['lon_v'] = lon_v
Gcut['lat_v'] = lat_v
Gcut['mask_v'] = mask_v
Gcut['lon_psi'] = lon_psi
Gcut['lat_psi'] = lat_psi

# then find the underlying plaid grid (x,y)
[x_rho,y_rho]=wfun.wcofs_lonlat_2_xy(lon_rho,lat_rho);
[x_u,y_u]=wfun.wcofs_lonlat_2_xy(lon_u,lat_u);
[x_v,y_v]=wfun.wcofs_lonlat_2_xy(lon_v,lat_v);
GG = dict()
GG['x_rho'] = x_rho
GG['y_rho'] = y_rho
GG['x_u'] = x_u
GG['y_u'] = y_u
GG['x_v'] = x_v
GG['y_v'] = y_v

# get dict of interpolants for each station
itp_dict = wfun.get_itp_dict(sta_dict, sta_list, GG, Gcut)

# creat lists of variables to get
v0_list = ['h', 'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v', 'angle']
v1_list = ['ocean_time']
v2_list = ['zeta']
v3_list_rho = ['salt', 'temp', 'u', 'v']
v3_list_w = []

# start by getting info from the first file
ds = nc.Dataset(fn_list[0])
dsg = nc.Dataset(grid_fn)

# generate dictionaries of long names and units
V_long_name = dict()
V_units = dict()
v_all_list = v1_list + v2_list + v3_list_rho + v3_list_w
for vv in v0_list:
    try:
        V_long_name[vv] = dsg.variables[vv].long_name
    except:
        V_long_name[vv] = ''
    try:
        V_units[vv] = dsg.variables[vv].units
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

# start netcdf files
for sta_name in sta_list:
    out_fn = out_fn_dict[sta_name]
    mfun.start_netcdf(out_fn, S['N'], NT, v0_list, v1_list, v2_list,
        v3_list_rho, v3_list_w, V_long_name, V_units)
        
# save static variables
for sta_name in sta_list:
    out_fn = out_fn_dict[sta_name]
    print('Initializing ' + out_fn)
    Xi0, Yi0, Xi1, Yi1, Aix, Aiy = itp_dict[sta_name]
    foo = nc.Dataset(out_fn, 'a')
    for vv in v0_list:
        xi01, yi01, aix, aiy = wfun.get_its(dsg, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
        vvtemp = Gcut[vv][yi01, xi01].squeeze()
        foo[vv][:] =   ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
    Lon = sta_dict[sta_name][0]
    Lat = sta_dict[sta_name][1]
    foo = nc.Dataset(out_fn)
    test_lon_rho = foo['lon_rho'][:]
    test_lat_rho = foo['lat_rho'][:]
    test_h = foo['h'][:]
    if False:
        # check on the interpolation, RESULT: it works as expected
        print(' Station Lon = %0.5f, Station Lat = %0.5f' % (Lon, Lat))
        print('    Test Lon = %0.5f,    Test Lat = %0.5f' % (test_lon_rho, test_lat_rho))
    foo.close()
ds.close()
dsg.close()
# END OF INITIALIZATION

# EXTRACT TIME-DEPENDENT FIELDS
count = 0
for fn in fn_list:
    ds = nc.Dataset(fn)
    print(' working on %d of %d' % (count, NT))
    sys.stdout.flush()
    for sta_name in sta_list:
        out_fn = out_fn_dict[sta_name]
        foo = nc.Dataset(out_fn, 'a')
        Xi0, Yi0, Xi1, Yi1, Aix, Aiy = itp_dict[sta_name]
        for vv in v1_list:
            vtemp = ds[vv][:].squeeze()
            foo[vv][count] = vtemp
        for vv in v2_list:
            xi01, yi01, aix, aiy = wfun.get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
            vvtemp = ds[vv][:, yi01, xi01].squeeze()
            vtemp = ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
            foo[vv][count] = vtemp
        for vv in v3_list_rho:
            xi01, yi01, aix, aiy = wfun.get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
            vvtemp = ds[vv][:, :, yi01, xi01].squeeze()
            vtemp = ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
            foo[vv][count,:] = vtemp
        for vv in v3_list_w:
            xi01, yi01, aix, aiy = wfun.get_its(ds, vv, Xi0, Yi0, Xi1, Yi1, Aix, Aiy)
            vvtemp = ds[vv][:, :, yi01, xi01].squeeze()
            vtemp = ( aiy*((aix*vvtemp).sum(-1)) ).sum(-1)
            foo[vv][count,:] = vtemp
        foo.close()
    count += 1

    ds.close()
# END OF EXTRACTING TIME-DEPENDENT FIELDS

# create z_rho and z_w (has to be done after we have zeta)
# also rotate velocity to be u = E-W, v = N-S
for sta_name in sta_list:
    out_fn = out_fn_dict[sta_name]
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

# PLOTTING
if False:
    import matplotlib.pyplot as plt
    plt.close('all')

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)

    hm = Gcut['h']
    hm[Gcut['mask_rho']==False] = np.nan
    cs = ax.pcolormesh(Gcut['lon_psi'], Gcut['lat_psi'], -hm[1:-1, 1:-1], cmap='rainbow', vmin=-1000, vmax=200)
    fig.colorbar(cs)

    for sta_name in sta_list:
        Lon = sta_dict[sta_name][0]
        Lat = sta_dict[sta_name][1]
        ax.plot(Lon,Lat,'ok')
        ax.text(Lon, Lat, sta_name)

    ax.grid(True)
    ax.set_title(Ldir['gtagex'])

    ax.axis([-126, -122, 44, 49])
    pfun.dar(ax)

    plt.show()


