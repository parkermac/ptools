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

do_movie = True

indir = os.path.abspath('../../ptools_output/bernoulli_average') + '/'
in_fn = indir + 'surface_lp.nc'

ds = nc.Dataset(in_fn)

# gather fields
xp = ds['lon_psi'][:]
yp = ds['lat_psi'][:]
xr = ds['lon_rho'][:]
yr = ds['lat_rho'][:]
ot = ds['ocean_time'][:]
nt = len(ot)

if True:
    aa = [-126, -122.1, 47, 51]
    xlist = [-126, -125, -124, -123]
    ylist = [47, 48, 49, 50, 51]
    vmin = -1
    vmax = 1
    aloc = 'upper right'
    outdir = os.path.abspath('../../ptools_output/bernoulli_movie_large') + '/'
else:    
    aa = [-124, -122.1, 47, 49]
    xlist = [-124, -123]
    ylist = [47, 48, 49]
    vmin = -.5
    vmax = .5
    aloc = 'lower left'
    outdir = os.path.abspath('../../ptools_output/bernoulli_movie_small') + '/'
Lfun.make_dir(outdir, clean=True)

# get indices for trimming
ix0 = zfun.find_nearest_ind(xr[0,:], aa[0])
ix1 = zfun.find_nearest_ind(xr[0,:], aa[1])
iy0 = zfun.find_nearest_ind(yr[:,0], aa[2])
iy1 = zfun.find_nearest_ind(yr[:,0], aa[3])

xxp = xp[iy0-1:iy1+1,ix0-1:ix1+1]
yyp = yp[iy0-1:iy1+1,ix0-1:ix1+1]

plt.close('all')
fs=16
plt.rc('font', size=fs)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
cmap = 'Spectral_r'

if do_movie:
    day_list = [item for item in range(len(ot))]
    # valid values are 0-89 so range(90) works
else:
    day_list = [0]
for day in day_list:
    zr = ds['zeta'][day,iy0:iy1+1,ix0:ix1+1].squeeze()
    uuvv = ds['uuvv'][day,iy0:iy1+1,ix0:ix1+1].squeeze()
    T = Lfun.modtime_to_datetime(ot[day])
    time_format = '%Y.%m.%d\n%H:%M'    
    Tstr = T.strftime(time_format)
    zr = zr - zr.mean()
    # really we should do an area-weighted mean but this is close enough

    # PLOTTING
    fig = plt.figure(figsize=(17,8))
    
    ax = fig.add_subplot(131)
    cs = ax.pcolormesh(xxp,yyp,9.8*zr, vmin=vmin, vmax=vmax, cmap=cmap)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xticks(xlist)
    ax.set_yticks(ylist)
    ax.set_title(r'$g\eta\ [m^{2}s^{-2}]$')
    ax.text(.05,.1,Tstr, transform=ax.transAxes)
    
    ax = fig.add_subplot(132)
    cs = ax.pcolormesh(xxp,yyp,0.5*uuvv, vmin=vmin, vmax=vmax, cmap=cmap)
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
    cs = ax.pcolormesh(xxp,yyp,9.8*zr + 0.5*uuvv, vmin=vmin, vmax=vmax, cmap=cmap)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xticks(xlist)
    ax.set_yticks([])
    ax.set_title(r'$g\eta+(u^{2} + v^{2})/2$')
    
    fig.tight_layout()
    if do_movie:
        nouts = ('0000' + str(day))[-4:]
        outname = 'plot_' + nouts + '.png'
        outfile = outdir + outname
        fig.savefig(outfile)
        plt.close()
    else:
        plt.show()
        
plt.rcdefaults()

ds.close()


if do_movie:    
    ff_str = ("ffmpeg -r 8 -i " + 
        outdir+"plot_%04d.png -vcodec libx264 -pix_fmt yuv420p -crf 25 "
        +outdir+"movie.mp4")
    os.system(ff_str)
