# -*- coding: utf-8 -*-
"""
Created on Mon May 23 18:13:32 2016

@author: PM5

Plotting a sphere for fun (and eventually for the tractive forces).
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

def get_tractive_scale():
    # distance from center of earth to center of moon [m]
    r = 384000e3
    # radius of the earth [m]
    r_e = 6371e3
    # mass of the moon [kg]
    M = 7.35e22
    # mass of the earth [kg]
    E = 5.972e24
    # acceleration of gravity [m s-2]
    g = 9.8
    # scale of the tractive force [m s-2]
    a = g*(M/E)*(r_e/r)**3
    return a

def get_TF(a, l,d):
    # returns vectors of the components of the tractive force
    # over the course of a lunar day
    # at latitude l [rad], for lunar declination d [rad]
    # longitude (as a way of telling time in lunar days)
    Z = np.linspace(0, 2*np.pi, 360)
    TFx = (3/2)*a*( np.sin(l)*np.sin(2*d)*np.sin(Z)
                   - np.cos(l)*np.cos(d)**2*np.sin(2*Z) )
    TFy = (3/2)*a*( -(1/2)*np.sin(2*l)*(1 - 3*np.sin(d)**2)
                   - np.cos(2*l)*np.sin(2*d)*np.cos(Z)
                   - (1/2)*np.sin(2*l)*np.cos(d)**2*np.cos(2*Z))
    return TFx, TFy


#%% plotting
plt.close()

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d')

r = 10
lat = np.linspace(-np.pi/2, np.pi/2, 100)
slice_rad = np.pi/6

lon = np.linspace(-np.pi, -slice_rad/2, 100)
x, y, z = get_xyz_mesh(lon, lat, r)
ax.plot_surface(x, y, z,
                rstride=4, cstride=4,
                color='y', linewidth=0, shade=True,
                alpha=.3)

lon = np.linspace(-slice_rad/2, slice_rad/2, 100)
x, y, z = get_xyz_mesh(lon, lat, r)
ax.plot_surface(x, y, z,
                rstride=4, cstride=4,
                color='r', linewidth=0, shade=True,
                alpha=.5)

lon = np.linspace(slice_rad/2, np.pi, 100)
x, y, z = get_xyz_mesh(lon, lat, r)
ax.plot_surface(x, y, z,
                rstride=4, cstride=4,
                color='y', linewidth=0, shade=True,
                alpha=.3)

lon = np.linspace(-np.pi, np.pi)
lat = 0*lon
X, Y, Z = get_xyz(lon, lat, r)
ax.plot(X, Y, Z, '-g', alpha=.4)

A = get_tractive_scale()
tf_dec_deg = 0
tf_dec = tf_dec_deg*np.pi/180
for tf_lat in np.pi*np.arange(-90,105,15)/180:
    TFx, TFy = get_TF(A, tf_lat, tf_dec)
    beta = np.pi/2 - tf_lat
    dy = (TFx/A)
    dx = -np.cos(beta) * (TFy/A)
    dz = np.sin(beta) * (TFy/A)
    x0, y0, z0 = get_xyz(0, tf_lat, r)
    X = x0 + dx
    Y = y0 + dy
    Z = z0 + dz
    ax.plot(X, Y, Z, '-k')
    if False:
        # Add a time marker
        # (makes lines less black for some reason)
        ii = 0
        ax.scatter(X[ii], Y[ii], Z[ii], s=5, c='k')

if False:
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

ax.set_title('Tractive Force for Lunar Declination = ' +
         str(tf_dec_deg) + '$^{o}$')

ax.azim = -18
ax.elev = 14

plt.show()

