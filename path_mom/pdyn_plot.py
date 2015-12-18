"""
Code to analyze the dynamics of flow along particle trajectories.

NOTE: we still need to do some interpolation in t, cs, lat, lon
before storing the results.
This may increase to runtime by a factor of 2**4 = 16.
Or we could just use the nearest of each
"""

# USER: set values
pmdir = '/Users/PM3/Documents/'
coast_file = pmdir + 'tools_data/geo_data/coast/pnw_coast_combined.mat'

# IMPORTS
import matplotlib.pyplot as plt
import numpy as np
import zfun; reload(zfun)

# make a list of roms files
rpath = pmdir + 'roms/output/D2005_dia/'
fn = rpath + 'ocean_dia_1800.nc'

# get the roms grid info
G, S = zfun.get_basic_info(fn, getT=False)
lon = G['lon_rho']
lat = G['lat_rho']
zbot = -G['h']

import cPickle
pin = open(pmdir + 'tools_output/pydev_out/p75.pkl', 'rb')
p = cPickle.load(pin)
pin.close()

px = p['lon']
py = p['lat']
NT, NP = px.shape

# determine distance from a point and then use that to decide
# winners and losers
# center point of circle (from the matlab code)
lon0 = -1.246791925614623e+02;
lat0 = 48.482319275124404;
clat = np.cos(np.deg2rad(lat0));
earth_rad = 6371e3 # m
x = earth_rad * clat * np.deg2rad(px - lon0)
y = earth_rad * np.deg2rad(py - lat0)
dist = np.sqrt(x**2 + y**2)

# also work with initial positions
x0 = x[0, :]
y0 = y[0, :]
z0 = p['z'][0, :]
zb0 = - p['H'][0, :]
deepMask = zb0 == zb0.min()
idm = np.where(deepMask)[0][0]
r0 = np.sqrt(x0[0]**2 + y0[0]**2) # radius in meters
a0 = np.arctan2(y0, x0) # radians
mNeg = a0 < 0.
mPos = a0 >= 0.
a0[mNeg] = np.pi + a0[mNeg]
a0[mPos] = a0[mPos] - np.pi
a00 = a0 # - a0[idm]
d00 = a00 * r0


# get the coastline
import scipy.io
cmat = scipy.io.loadmat(coast_file)
cstlon = cmat['lon'].squeeze()
cstlat = cmat['lat'].squeeze()

plt.close()

fig = plt.figure(figsize=(16,8))
# map
ax = fig.add_subplot(121)
plt.rcParams['contour.negative_linestyle'] = 'solid'
ax.contour(lon, lat, zbot, [-2000, -1000, -500, -200, -100], colors='g')
ax.axis([-126.5, -123.5, 47, 49.5])
#ax.set_xlim(lon.min(), lon.max())
#ax.set_ylim(lat.min(), lat.max())
ax.plot(lon0, lat0,'oy', markersize=10)
iwin = []
for ii in range(NP):
    if np.nanmin(dist[:, ii]) < 10e3:
        iend = np.where(dist[:, ii] < 10e3)[0][0]
        if p['z'][iend, ii] < -70.0:
            ax.plot(px[:iend+1, ii], py[:iend+1, ii], '-r', linewidth=1)
            ax.plot(px[0, ii], py[0, ii], '*r', markersize=15)
            iwin.append(ii)      
    else:
        ax.plot(px[:, ii], py[:, ii], '-k', linewidth=.2)
        ax.plot(px[0, ii], py[0, ii], 'oc', markersize=3)
     
zfun.dar(ax)
# add coastline
ax.plot(cstlon, cstlat, '-k', linewidth=.5)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Bathymetry and Tracks')

ax = fig.add_subplot(122)
cs0 = p['cs'][0, :] # all initial cs values
css = set(np.round(cs0*100)) # get just unique values
# need to round because values are not exact
nloc = len(css) # number of unique initial cs values
ax.plot(d00[:NP/nloc]/1000, zb0[:NP/nloc], '-', color='0.5', linewidth=3)
ax.plot(d00/1000, z0, 'oc', markersize=3)
ax.plot(d00[iwin]/1000, z0[iwin], '*r', markersize=15)

plt.show()

