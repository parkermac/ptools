"""
Code to explore visualization of 2D currents.
"""

# setup
import os; import sys
alp = os.path.abspath('/Users/PM3/Documents/LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart(alp)
import zfun; reload(zfun)
import matfun; reload(matfun)

import numpy as np
import matplotlib.pyplot as plt

dir0 = '/Users/PM3/Documents/roms/output/D2005_his/'
nhis = 1800
fn = dir0 + 'ocean_his_' + str(nhis) + '.nc'

G, S, T = zfun.get_basic_info(fn)

# coastline file
coast_file = Ldir['data'] + 'coast/pnw_coast_combined.mat'

# GET DATA
G, S, T = zfun.get_basic_info(fn)
import netCDF4 as nc   
ds = nc.Dataset(fn,'r')
h = G['h']
salt = ds.variables['salt'][0, -1, 1:-1, 1:-1].squeeze()
temp = ds.variables['temp'][0, -1, 1:-1, 1:-1].squeeze()

lon = G['lon_psi']
lat = G['lat_psi']

# get the coastline
cmat = matfun.loadmat(coast_file)

# plotting
plt.close()
fig = plt.figure(figsize=(14, 8))

# surface salinity
ax = fig.add_subplot(121)
cmap = plt.get_cmap(name='jet')
cs = ax.pcolormesh(lon, lat, salt, vmin=28, vmax=34,  cmap = cmap)
ax.contour(G['lon_rho'], G['lat_rho'], h,
    [100, 200, 500, 1000, 2000, 3000], colors='g')
ax.axis([lon.min(), lon.max(), lat.min(), lat.max()])
# add coastline
ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
zfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Surface Salinity')
fig.colorbar(cs)

# put info on plot
xt, yt = zfun.labxy(ax, 'lr')
f_string = T['tm'].strftime('%Y.%m.%d')
ax.text(xt, yt + .6, f_string, horizontalalignment='right')
ax.text(xt, yt + .4, T['tm'].strftime('%Y-%m-%d'), horizontalalignment='right')
ax.text(xt, yt + .2, T['tm'].strftime('%H:%M'), horizontalalignment='right')
ax.text(xt, yt, 'Save = ' + str(nhis), horizontalalignment='right')

# surface temperature
ax = fig.add_subplot(122)
cmap = plt.get_cmap(name='rainbow')
cs = ax.pcolormesh(lon, lat, temp, vmin=8, vmax=14,  cmap = cmap)
ax.contour(G['lon_rho'], G['lat_rho'], h,
    [100, 200, 500, 1000, 2000, 3000], colors='g')
ax.axis([lon.min(), lon.max(), lat.min(), lat.max()])
# add coastline
ax.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
zfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Surface Temperature')
fig.colorbar(cs)

# add velocity vectors

# get velocity
u = ds.variables['u'][0, -1, :, :].squeeze()
v = ds.variables['v'][0, -1, :, :].squeeze()

# set masked values to 0
um = u.data; um[G['mask_u']==False] = 0
vm = v.data; vm[G['mask_v']==False] = 0

# full coordinate arrays
Xu = G['lon_u']
Yu = G['lat_u']
Xv = G['lon_v']
Yv = G['lat_v']

# axis vectors (plaid grids)
xu = Xu[0,:]
yu = Yu[:,0]
xv = Xv[0,:]
yv = Yv[:,0]

def step_ahead(x, y, u, v, dt):
    # all inputs except dt neet to be ndarrays of the same size
    # x and y are degrees lon and lat
    RE = 6371.0e3
    clat = np.cos(np.deg2rad(y))
    x1 = x + np.rad2deg(u*dt/(RE*clat))
    y1 = y + np.rad2deg(v*dt/(RE))   
    return x1, y1

x1, y1 = step_ahead(Xu, Yv, um, vm, 3600.)

# you can use np.take(array, indices of flattened array) to get selected
# elements

# what I want to build is a function
# interp2d_plaid_to_scatter which returns f(xi, yi) for arbitrary positions
# xi and yi (each can be a vector) starting from f(x,y) where the grid is plaid.




plt.show()


