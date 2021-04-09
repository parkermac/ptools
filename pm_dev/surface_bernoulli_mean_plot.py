"""
Exploring the tidally averaged surface Bernoulli function
from a LiveOcean layer extraction.
"""

# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from time import time

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
import zfun

sys.path.append(os.path.abspath('../../LiveOcean/plotting'))
import pfun

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

Ldir = Lfun.Lstart()

do_all = True

indir = os.path.abspath('../../ptools_output/bernoulli_average') + '/'
in_fn = indir + 'surface_lp.nc'

ds = nc.Dataset(in_fn)

# gather fields
xp = ds['lon_psi'][:]
yp = ds['lat_psi'][:]
xr = ds['lon_rho'][:]
yr = ds['lat_rho'][:]
ot = ds['ocean_time'][:]

if True:
    aa = [-126, -122.1, 47, 51]
    xlist = [-126, -125, -124, -123]
    ylist = [47, 48, 49, 50, 51]
    vmin = -1
    vmax = 1
    aloc = 'upper right'
    fn = 'bernoulli_mean_large.png'
else:    
    aa = [-124, -122.1, 47, 49]
    xlist = [-124, -123]
    ylist = [47, 48, 49]
    vmin = -.5
    vmax = .5
    aloc = 'lower left'
    fn = 'bernoulli_mean_small.png'
out_fn = indir + fn

# get indices for trimming
ix0 = zfun.find_nearest_ind(xr[0,:], aa[0])
ix1 = zfun.find_nearest_ind(xr[0,:], aa[1])
iy0 = zfun.find_nearest_ind(yr[:,0], aa[2])
iy1 = zfun.find_nearest_ind(yr[:,0], aa[3])

xxp = xp[iy0-1:iy1+1,ix0-1:ix1+1]
yyp = yp[iy0-1:iy1+1,ix0-1:ix1+1]
xxr = xr[iy0-1:iy1+1,ix0-1:ix1+1]
yyr = yr[iy0-1:iy1+1,ix0-1:ix1+1]

if do_all:
    pass
else:
    ot = ot[:3]
nt = len(ot)
    
day_list = [item for item in range(len(ot))]
for day in day_list:
    zr = ds['zeta'][day,iy0:iy1+1,ix0:ix1+1].squeeze()
    uuvv = ds['uuvv'][day,iy0:iy1+1,ix0:ix1+1].squeeze()
    u = ds['u'][day,iy0:iy1+1,ix0:ix1+1].squeeze()
    v = ds['v'][day,iy0:iy1+1,ix0:ix1+1].squeeze()
    zr = zr - zr.mean()
    if day == 0:
        zrm = zr/nt
        uuvvm = uuvv/nt
        um = u/nt
        vm = v/nt
    else:
        zrm += zr/nt
        uuvvm += uuvv/nt
        um += u/nt
        vm += v/nt
        
ds.close()

time_format = '%Y.%m.%d'
T0 = Lfun.modtime_to_datetime(ot[0])
Tstr0 = T0.strftime(time_format)
T1 = Lfun.modtime_to_datetime(ot[-1])
Tstr1 = T1.strftime(time_format)

# PLOTTING
plt.close('all')
fs=16
plt.rc('font', size=fs)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
cmap = 'Spectral_r'

if False:
    # Bernoulli plot
    fig = plt.figure(figsize=(17,8))

    ax = fig.add_subplot(131)
    cs = ax.pcolormesh(xxp,yyp,9.8*zrm, vmin=vmin, vmax=vmax, cmap=cmap)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xticks(xlist)
    ax.set_yticks(ylist)
    ax.set_title(r'$g\eta\ [m^{2}s^{-2}]$')
    ax.text(.05,.1,Tstr0, transform=ax.transAxes)
    ax.text(.05,.05,Tstr1, transform=ax.transAxes)

    ax = fig.add_subplot(132)
    cs = ax.pcolormesh(xxp,yyp,0.5*uuvvm, vmin=vmin, vmax=vmax, cmap=cmap)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xticks(xlist)
    ax.set_yticks([])
    ax.set_title(r'$(u^{2} + v^{2})/2$')
    # Inset colorbar
    cbaxes = inset_axes(ax, width="4%", height="40%", loc=aloc, borderpad=3) 
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')

    ax = fig.add_subplot(133)
    cs = ax.pcolormesh(xxp,yyp,9.8*zrm + 0.5*uuvvm, vmin=vmin, vmax=vmax, cmap=cmap)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xticks(xlist)
    ax.set_yticks([])
    ax.set_title(r'$g\eta+(u^{2} + v^{2})/2$')
    
else:
    # speed plot
    fig = plt.figure(figsize=(8,12))
    
    ax = fig.add_subplot(111)
    cs = ax.pcolormesh(xxp,yyp,np.sqrt(um*um + vm*vm), vmin=0, vmax=0.6, cmap=cmap)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xticks(xlist)
    ax.set_yticks(ylist)
    ax.set_title(r'Speed of Mean Currents $[m\ s^{-1}]$')
    nskp = 5
    ax.quiver(xxr[::nskp,::nskp], yyr[::nskp,::nskp], um[::nskp,::nskp], vm[::nskp,::nskp],
        units='y', scale=10, scale_units='y', color='gray', pivot='middle')
    # Inset colorbar
    cbaxes = inset_axes(ax, width="4%", height="40%", loc=aloc, borderpad=4) 
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    

fig.tight_layout()
plt.show()
plt.savefig(out_fn)
plt.rcdefaults()


