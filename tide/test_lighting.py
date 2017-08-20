"""
Code to test lighting of surfaces.

http://physicalmodelingwithpython.blogspot.com/2015/08/illuminating-surface-plots.html
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D     # Import 3D plotting tools.
from scipy.special import jn                # Import Bessel function.

def get_xyz_mesh(lon, lat, r):
    # makes a Cartesian mesh on a sphere
    # lon, lat are numpy vectors
    # r is a scalar (the radius)
    x = r * np.outer(np.cos(lon), np.cos(lat))
    y = r * np.outer(np.sin(lon), np.cos(lat))
    z = r * np.outer(np.ones(np.size(lon)), np.sin(lat))
    return x,y,z
    
# # Define grid of points.
# points = np.linspace(-10, 10, 51)
# X, Y = np.meshgrid(points, points)
# R = np.sqrt(X**2 + Y**2)
# Z = jn(0,R)

r = 10
lat = np.linspace(-np.pi/2, np.pi/2, 100)
# first part of the background sphere
lon = np.linspace(-np.pi, np.pi, 100)
X,Y,Z = get_xyz_mesh(lon, lat, r)
# ax.plot_surface(x, y, z,
#                 rstride=4, cstride=4,
#                 color='b', linewidth=0, shade=True,
#                 alpha=.4)

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

green = np.array([0,1.0,0])
green_surface = light.shade_rgb(rgb * green, Z)

ax = Axes3D(plt.figure())
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0, antialiased=False,
                facecolors=green_surface)
                
plt.show()

