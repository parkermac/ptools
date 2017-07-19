#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 15 15:58:43 2017

@author: PM5

Numerically simulate frictional tides in a channel of irregular
cross-section.

"""

import numpy as np
import matplotlib.pyplot as plt

# channel parameters
L = 600e3
x = np.linspace(-L, 0, 200)
xx = np.linspace(-2, 2, 200)
ovec = np.ones_like(x)
H = 150*ovec - 100 * np.exp(-xx**2)

# other parameters
g = 9.8
Th = 12.42
om = 2*np.pi/(Th*3600)
c = np.sqrt(g*H)
Cd = 5e-3
U = 1
r = Cd * U / H

# calculate complex wavenember K(x)
K = (om/c) * np.sqrt(1 + (r/om)*1j)

# specify complex amplitude of ocean tide
A = 1
B = 0

# calculate complex amplitude of response (both functions of x)
a = ((A + B*1j) / (np.exp(1j*(-K*L)) + np.exp(1j*(K*L)))).real
b = ((A + B*1j) / (np.exp(1j*(-K*L)) + np.exp(1j*(K*L)))).imag

# function to calculate eta(t)
def get_eta(a,b,K,x,om,t):
    eta = ( (a + b*1j) * ( np.exp(1j*(K*x-om*t)) + np.exp(1j*(-K*x-om*t)) ) ).real
    return eta

# plotting
xkm = x/1e3
Lkm = L/1e3
plt.close('all')
fig = plt.figure(figsize=(13,8))

ax = fig.add_subplot(2,1,1)
for t in np.linspace(0, 86400*(Th/12), 2*12 + 1):    
    ax.plot(xkm, get_eta(a,b,K,x,om,t))
ax.set_xlim(-Lkm, 0)
ax.set_ylabel('eta (m)')
    
ax = fig.add_subplot(2,1,2)
ax.plot(xkm, -H)
ax.set_xlim(-Lkm, 0)
ax.set_ylim(top=0)
ax.set_xlabel('x (km)')
ax.set_ylabel('bathymetry z (m)')

plt.show()