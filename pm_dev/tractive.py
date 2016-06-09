# -*- coding: utf-8 -*-
"""
Created on Fri May  6 16:56:56 2016

@author: PM5

Code to calculate the lunar tractive force.
"""

import numpy as np
import matplotlib.pyplot as plt

# latitude
lat_deg_vec = np.array([85, 60, 30, 0])
# declination of the moon
dec_deg_vec = np.array([0, 15, 28])

# longitude (as a way of telling time in lunar days)
Z = np.linspace(0, 2*np.pi, 360)
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
# c is a latitude-like angle starting from zero at the
# axis pointing to the moon [rad]
c = np.linspace(0, np.pi,1000)
# components of the tractive force (relative to c)
TFv = a*(2*np.cos(c)**2 - np.sin(c)**2) # vertical (dynamically unimportant)
TFh = (3/2)*a*np.sin(2*c) # horizontal
TF = (3/2)*a

def get_TF(a,l,d):
    # returns vectors of the components of the tractive force
    # over the course of a lunar day
    # at latitude l [rad], for lunar declination d [rad]
    TFx = (3/2)*a*( np.sin(l)*np.sin(2*d)*np.sin(Z)
                   - np.cos(l)*np.cos(d)**2*np.sin(2*Z) )
    TFy = (3/2)*a*( -(1/2)*np.sin(2*l)*(1 - 3*np.sin(d)**2)
                   - np.cos(2*l)*np.sin(2*d)*np.cos(Z)
                   - (1/2)*np.sin(2*l)*np.cos(d)**2*np.cos(2*Z))
    return TFx, TFy

#plotting

plt.close('all')

# the tractive force (and the vertical component)
# relative to the axis pointing toward the moon
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
ax.plot(180*c/np.pi, TFh/TF, '-k')
ax.plot(180*c/np.pi, TFv/TF, '--k')
ax.legend(('Horizontal','Vertical'))
ax.set_title('max(TF) = %.3g [m s-2]' % (np.max(TF)))
ax.set_xlabel('c (degrees from axis pointing to Moon)')
ax.set_ylabel('TF/max(TF)')
ax.set_xlim(0, 180)
ax.grid()

# the horizontal tractive force at various latitudes and lunar declinations
# over the course of 24 lunar hours
NR = len(lat_deg_vec)
NC = len(dec_deg_vec)
fig, axes = plt.subplots(nrows=NR, ncols=NC,
                         figsize=(3*NC, 3*NR), squeeze=False)
cc = 0
for l_deg in lat_deg_vec:
    for d_deg in dec_deg_vec:
        l = l_deg*np.pi/180
        d = d_deg*np.pi/180
        TFx, TFy = get_TF(a, l, d)
        ir = int(np.floor(cc/NC))
        ic = int(cc - NC*ir)
        ax = axes[ir, ic]
        ax.plot(TFx/TF, TFy/TF)
        if ir == 0:
            ax.set_title('Lunar Dec. = %d' %
                        (d_deg))
        if ic == 0:
            ax.text(-.9, .85, 'Latitude = ' + str(l_deg))
            ax.set_ylabel('TFy/max(TF)')
        if ir == NR-1:
            ax.set_xlabel('TFx/max(TF)')
        ax.set_xlim(-1.1, 1.1)
        ax.set_ylim(-1.1, 1.1)
        for ii in range(24):
            ax.text(TFx[ii*15]/TF, TFy[ii*15]/TF - ii/200, str(ii))
        ax.set_aspect('equal')
        ax.grid()
        cc += 1

plt.show()

