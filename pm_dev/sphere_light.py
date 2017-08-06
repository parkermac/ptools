# -*- coding: utf-8 -*-
"""
Plotting a 3D sphere with lighting

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def get_xyz_mesh(lon, lat, r):
    # makes a Cartesian mesh on a sphere
    # lon, lat are numpy vectors
    # r is a scalar (the radius)
    x = r * np.outer(np.cos(lon), np.cos(lat))
    y = r * np.outer(np.sin(lon), np.cos(lat))
    z = r * np.outer(np.ones(np.size(lon)), np.sin(lat))
    return x,y,z

def get_xyz(lon, lat, r):
    # makes x,y,z vectors on a sphere
    # lon, lat are numpy vectors (equal length)
    # r is a scalar (the radius)
    x = r * np.cos(lon) * np.cos(lat)
    y = r * np.sin(lon) * np.cos(lat)
    z = r * np.sin(lat)
    return x,y,z

#%% plotting
plt.close()

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection='3d')

r = 10
lat = np.linspace(-np.pi/2, np.pi/2, 100)
slice_rad = np.pi/6

lon = np.linspace(-np.pi, np.pi, 100)
x, y, z = get_xyz_mesh(lon, lat, r)

from matplotlib.colors import LightSource
from matplotlib import cm
light = LightSource(45,45)
#illuminated_surface = light.shade(z, cmap=cm.rainbow)

rgb = np.ones((z.shape[0], z.shape[1], 3))
illuminated_surface = light.shade_rgb(rgb, z)
#ax = Axes3D(plt.figure())
ax.plot_surface(x,y,z, rstride=4, cstride=4, linewidth=0, antialiased=False, shade=True)
                #facecolors=illuminated_surface)

# ax.plot_surface(x, y, z,
#                 rstride=4, cstride=4,
#                 color='y', linewidth=0, shade=True,
#                 alpha=1, facecolors=illuminated_surface)

lon = np.linspace(-np.pi, np.pi)
lat = 0*lon
X, Y, Z = get_xyz(lon, lat, r)



# ax.plot(X, Y, Z, '-g', alpha=.4)


if True:
    # for development
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    scl = 1.2
else:
    # much cleaner display
    scl = .8
    ax.set_axis_off()

ax.set_xlim(-scl*r, scl*r)
ax.set_ylim(-scl*r, scl*r)
ax.set_zlim(-scl*r, scl*r)
ax.set_aspect('equal')

ax.azim = -18
ax.elev = 14

plt.show()

