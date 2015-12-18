"""
3D plot of the volume of Puget Sound.
"""

# setup
import os; import sys
alp = os.path.abspath('/Users/PM3/Documents/LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import zfun; reload(zfun)
import matfun; reload(matfun)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

dir0 = '/Users/PM3/Documents/tools_data/geo_data/topo/psdem/'

fn = dir0 + 'psdem_2005_plaid_183m.nc'

import netCDF4 as nc

ds = nc.Dataset(fn, 'r')

if False:
    for vn in ds.variables:
        print ds.variables[vn]
    
z = ds.variables['elev'][:]
lon = ds.variables['lon'][:]
lat = ds.variables['lat'][:]

RE = 6371.0e3
clat = np.cos(np.deg2rad(lat.mean()))
x = RE * np.deg2rad(lon - lon[0]) * clat
y = RE * np.deg2rad(lat - lat[0])
X, Y = np.meshgrid(x/1000,y/1000)

# convert z to a masked array
z[z>0] = np.nan
Z = np.ma.masked_invalid(z)
D = -Z


plt.close()

fig = plt.figure(figsize=(14,8))
ax = fig.add_subplot(111, projection='3d')

#ax.plot_surface(X, Y, D, rstride=5, cstride=5,
#    linewidth=0, cmap=cm.rainbow, vmin=0, vmax = 100)
    
ax.contour(X, Y, D, levels=np.arange(0,300,10), rstride=1, cstride=1,
    cmap=cm.rainbow, vmin=0, vmax = 100)

        
plt.show()




