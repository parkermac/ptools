"""
This code is to test out some ways to make distance and grounding calculations
faster.  Done as an example for Marissa Leatherman for the Copper project, and
to help me thing through the issues.

Note that the distance is along the track, so it includes all the tidal
back-and-forth.  It will be much longer than the as-the-crow-flies distance.

Code should be self-contained.

"""

import netCDF4 as nc
import numpy as np
import seawater as sw # (get from conda install seawater)
import matplotlib.pyplot as plt

fn = ('/Users/pm8/Documents/LiveOcean_output/tracks/'
        + 'leathermanPuy_3d_nosink/release_2018.01.01.nc')
        
ds = nc.Dataset(fn)
# variables are
# ['lon', 'lat', 'cs', 'ot', 'z', 'salt', 'temp', 'zeta', 'h', 'u', 'v', 'w']
# packed as (time, particle)

# get time in days since the release
ts = ds['ot'][:].data
td = (ts - ts[0])/86400

# goal is to find out how far each particle goes along its track before
# it sticks (defined as cs <= -0.99)

# calculate all distances
lon = ds['lon'][:].data
lat = ds['lat'][:].data
# we take .data because these are masked arrays
NT, NP = lon.shape
# pre-allocate a distance array
dist = np.zeros((NT, NP))
# dist is distance between points [km]
for p in range(NP):
    # add one particle at a time
    dist[1:,p], junk = sw.dist(lat[:,p], lon[:,p])
    # note that the distance arrays are one less in length than the positions
    # because they are between points
cdist = np.cumsum(dist, axis=0)
# cdist is the alongtrack distance from the start
# Note that the array returned by cumsum is the same size as the input array,
# and by the way we made "dist" the first point will have cdist = 0, as expected.

# finally set each track to be nan's after the particle has "stuck"
cs = ds['cs'][:].data
stuck = cs <= -.99 # a Boolean array that is True everyplace the stuck condition applies
not_stuck = ~stuck
istuck = np.argmin(not_stuck, axis = 0)
# argmin applied to a Boolean array gives the index of the first occurence of False,
# in this case it is for each column because I used axis = 0.
# by using not_stuck it means that the array argmin is searching has False where a
# particle should be stuck.

# make and array that has nan's after the particle gets stuck.
cdist_alt = cdist.copy()
for p in range(NP):
    cdist_alt[istuck[p]:,p] = np.nan

# plotting
plt.close('all')
fs = 14
plt.rc('font', size=fs)
fig = plt.figure(figsize=(8,8))

ax = fig.add_subplot(111)
for p in range(100):
    ax.plot(td,cdist[:,p], lw = 1, alpha = .5) # full path
    ax.plot(td,cdist_alt[:,p], lw = 2) # ends when it gets stuck
ax.set_xlabel('Time [days]')
ax.set_ylabel('Distance [km]')

plt.show()
plt.rcdefaults()

