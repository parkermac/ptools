"""
Code to test lighting of surfaces.

http://physicalmodelingwithpython.blogspot.com/2015/08/illuminating-surface-plots.html
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D     # Import 3D plotting tools.
from scipy.special import jn                # Import Bessel function.

#%% setup
import os
import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart(gridname='cascadia1', tag='base')
Ldir['gtagex'] = Ldir['gtag'] + '_' + 'lobio1'
import netCDF4 as nc
import zrfun
    
# # Define grid of points.
# points = np.linspace(-10, 10, 51)
# X, Y = np.meshgrid(points, points)
# R = np.sqrt(X**2 + Y**2)
# Z = jn(0,R)

indir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f2017.08.07/'
fn = 'ocean_his_0002.nc'
in_fn = indir + fn

G = zrfun.get_basic_info(in_fn, only_G=True)

X = G['lon_rho']
Y = G['lat_rho']
ds = nc.Dataset(in_fn)
ZZ = ds['zeta'][:].squeeze()
Z = ZZ.data
Z[ZZ.mask==True] = np.nan
# Get lighting object for shading surface plots.
from matplotlib.colors import LightSource

# Get colormaps to use with lighting object.
from matplotlib import cm

# Create an instance of a LightSource and use it to illuminate the surface.
light = LightSource(90, 45)
illuminated_surface = light.shade(Z, cmap=cm.coolwarm)

plt.close('all')

# Create 3D surface plot.

# Create a shaded green surface.

rgb = np.ones((Z.shape[0], Z.shape[1], 3))

clr = np.array([1.0, 0.7, 0.5]) # copper
clr_surface = light.shade_rgb(rgb * clr, Z)

fig = plt.figure(figsize=(14,8))
ax = Axes3D(fig)
ax.plot_surface(X, Y, 0.1*Z, rstride=2, cstride=2, linewidth=0, antialiased=False,
                facecolors=clr_surface)

yl = ax.get_ylim()
yav = (yl[0] + yl[1])/2
#ax.set_aspect(1/np.sin(np.pi*yav/180))
ax.set_aspect('equal')

#ax.set_zlim(-10, 10)

ax.azim = 165
ax.elev = 15

#ax.set_axis_off()


plt.show()

