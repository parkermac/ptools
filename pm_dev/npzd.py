#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:08:29 2017

@author: PM5

Code for time-dependent NPZD simulation.
Based on Banas etal. (2009) JGR, RISE

"""

import numpy as np
import matplotlib.pyplot as plt

def mu_i(E, N):
    # instantaneous phytoplankton growth rate
    mu0 = 2.2 # maximum instantaneous growth rate: d-1
    ks = 4.6 # half-satuaration for nitrate uptake: uM N
    alpha = 0.07 # initial slope of growth-ligh curve: (W m-2)-1 d-1    
    mu_i = mu0 * (N/(ks + N)) * (alpha * E)/np.sqrt(mu0**2 + alpha**2*E**2)
    return mu_i

def I(P):
    # zooplankton ingestion
    I0 = 4.8 # maximum ingestion rate: d-1
    Ks = 3 # half-satuaration for ingestion: uM N
    I = I0 * P**2 / (Ks**2 + P**2)
    return I

# all time will be in days
#
T = 20 # total time: d
dt = .01 # time step: d
NT = int(T/dt) # number of time steps

t = np.nan * np.ones(NT) # time d

N = np.nan * np.ones(NT) # nitrate uM N
P = np.nan * np.ones(NT) # phytoplankton uM N
Z = np.nan * np.ones(NT) # zooplankton uM N
D = np.nan * np.ones(NT) # detritus uM N

def Ef(t):
    E0 = 200 # max light: W m-2
    E = (E0/2) * (1 + np.cos(2*np.pi*t))
    return E

m = 0.1 # nongrazing mortality: d-1
eps = 0.3 # growth efficiency
xi = 2.0 # mortality: d-1 (uM N)-1
f_e = 0.5 # fraction of losses egested
r = 0.1 # remineralization rate d-1

N[0] = 30
P[0] = 1
Z[0] = .1
D[0] = 0
t[0] = 0
ii = 0

while ii < NT-1:    
    ii += 1    
    t[ii] = t[ii-1] + dt
    
    E = Ef(t[ii])
    
    Np = N[ii-1]    
    Pp = P[ii-1]
    Zp = Z[ii-1]
    Dp = D[ii-1]
    
    N[ii] = Np + dt*( -mu_i(E, Np)*Pp + (1-eps)*(1-f_e)*I(Pp)*Zp + r*Dp )
    P[ii] = Pp + dt*( mu_i(E, Np)*Pp - I(Pp)*Zp - m*Pp )
    Z[ii] = Zp + dt*( eps*I(Pp)*Zp - xi*Zp**2 )
    D[ii] = Dp + dt*( (1-eps)*f_e*I(Pp)*Zp + m*Pp + xi*Zp**2 - r*Dp)
    
# plotting

plt.close('all')
fig = plt.figure(figsize=(13,8))
ax = fig.add_subplot(1,1,1)

lw = 4
ax.plot(t, N, '-', color='royalblue', linewidth=lw)
ax.plot(t, P, '-', color='lightgreen', linewidth=lw)
ax.plot(t, Z, '-', color='pink', linewidth=lw)
ax.plot(t, D, '-', color='tan', linewidth=lw)
ax.legend(['Nitrate', 'Phytoplankton', 'Zooplankton', 'Detritus'])
ax.set_xlim(t[0], t[-1])
ax.grid()

ax.set_xlabel('Days')
ax.set_ylabel('$(\mu mol\ N\ L^{-1})$')

plt.show()

