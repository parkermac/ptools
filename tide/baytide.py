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
L = 550e3
nx = 100
x = np.linspace(-L, 0, nx)
xx = np.linspace(-2, 3, nx)
#ovec = np.ones(x)
H = np.tanh(xx) * (144 - 106)/2 + (144 + 106)/2
H = H  - 60 * np.exp(-xx**2)

# note: for H=50, and omega for M2, the wavelength is ~1000 km.

# other parameters
g = 9.8
Th = 12.42
om = 2*np.pi/(Th*3600)
c = np.sqrt(g*H)
Cd = 5e-3
U = 1
r = Cd * U / H

# calculate complex wavenember K(x)
def get_K(om, c, r):
    K = (om/c) * np.sqrt(1 + (r/om)*1j)
    return K
K = get_K(om, c, r)

# specify complex amplitude of ocean tide
gam = 1 + 0*1j

# find solution using matrix inversion
def get_alpha(K, L, gam):
    EL = np.exp(1j*K[0]*L)
    A = np.matrix([[1, -1], [1/EL, EL]])
    Y = np.matrix([0, gam]).T
    AI = np.linalg.inv(A)
    X = AI*Y
    #Xalt = np.linalg.solve(AA, Y) # identical to X   
    alpha_p = X[0,0]
    alpha_m = X[1,0]
    return alpha_p, alpha_m
alpha_p, alpha_m = get_alpha(K, L, gam)

# function to calculate eta(t)
def get_eta(alpha_p, alpha_m, K, x, om, t):
    eta = ( alpha_p*np.exp(1j*(K*x-om*t)) + alpha_p*np.exp(1j*(-K*x-om*t)) ).real
    return eta
    
# function to calculate amplitude and phase as a function of x
def get_amp_ph(alpha_p, alpha_m, K, om, x):
    KR = K.real
    KI = K.imag
    f_p = alpha_p * np.exp(-KI*x) * np.exp(1j*KR*x)
    f_m = alpha_m * np.exp(KI*x) * np.exp(-1j*KR*x)
    f = f_p + f_m
    a = f.real
    b = f.imag
    amp = np.sqrt(a**2 + b**2)
    ph = np.arctan2(b,a)
    return amp, ph, a, b

amp, ph, a, b = get_amp_ph(alpha_p, alpha_m, K, om, x)

# plotting

xkm = x/1e3
Lkm = L/1e3
plt.close('all')
fig = plt.figure(figsize=(13,8))

ax = fig.add_subplot(2,2,1)
for t in np.linspace(0, 86400*(Th/12), 2*12 + 1):    
    ax.plot(xkm, get_eta(alpha_p, alpha_m, K, x, om, t).T)
ax.set_xlim(-Lkm, 0)
ax.set_ylabel('eta (m)')
ax.set_title('Eta over a tidal cycle')
    
ax = fig.add_subplot(2,2,2)
# apmplitude
ax.plot(xkm, amp, '-r')
ax.set_xlim(-Lkm, 0)
ax.set_ylabel('Amplitude (m)', color='r')
ax.tick_params('y', colors='r')
# phase
ax2 = ax.twinx()
ax2.plot(xkm, ph/np.pi, '-b', label='Phase Lag / 180 deg')
ax2.set_xlim(-Lkm, 0)
ax2.set_ylim(0, 1)
ax2.set_ylabel('Phase Lag / 180 deg', color='b')
ax2.tick_params('y', colors='b')
ax2.grid()

# Bathymetry
ax3 = fig.add_subplot(2,2,3)
ax3.plot(xkm, -H)
ax3.set_xlim(-Lkm, 0)
ax3.set_ylim(top=0)
ax3.set_ylabel('Bathymetry (m)')
ax3.set_xlabel('x (km)')


# plot versus frequency

Th0 = 12.42
om0 = 2*np.pi/(Th0*3600)

om_vec = np.linspace(om0/100, 1.5*om0, 100)
amp_vec = np.nan * om_vec
ph_vec = np.nan*om_vec
ii = 0
for om1 in om_vec:
    K = get_K(om1, c, r)
    alpha_p, alpha_m = get_alpha(K, L, gam)
    amp, ph, a, b = get_amp_ph(alpha_p, alpha_m, K, om, x)
    amp_vec[ii] = amp[-1]
    ph_vec[ii] = ph[-1]
    ii += 1

mask = ph_vec < 0
ph_vec[mask] = ph_vec[mask] + 2*np.pi

ax = fig.add_subplot(2,2,4)
# apmplitude
ax.plot(om_vec/om0, amp_vec, '-r')
ax.set_xlim(0, om_vec[-1]/om0)
ax.set_ylabel('Amplitude (m)', color='r')
ax.tick_params('y', colors='r')
ax.grid()
# phase
ax2 = ax.twinx()
ax2.plot(om_vec/om0, ph_vec/np.pi, '-b')
ax2.set_xlim(0, om_vec[-1]/om0)
ax2.set_ylim(bottom=0)
ax2.set_ylabel('Phase Lag / 180 deg', color='b')
ax2.tick_params('y', colors='b')
ax.set_xlabel('Frequency / M2 Frequency')

# add observational data from Sutherland etal (2005)
# using station 5 for the ocean and Little River for the head of SoG
Cons = {'M2':12.42, 'S2':12, 'N2':12.66, 'K2':11.97,
        'K1':23.93, 'O1':25.82, 'P1':24.07, 'Q1':26.87}
Amps = [(94.4, 99.36), (25.9, 25.02), (18.6, 21.64), (5.9, 6.8),
        (43.3, 90.19), (25.6, 49.26), (13.2, 28.62), (5, 8.38)]
Phases = [(233.3, 32.87), (261.2, 61.6), (207.6, 5.42), (247.8, 62.56),
          (238, 287.03), (223.6, 263.94), (239.3, 285.67), (222.6, 257.2)]

ii = 0
for Con in Cons.keys():
    Om = 2*np.pi/(3600*Cons[Con])
    Amp = Amps[ii]
    dAmp = Amp[1]/Amp[0]
    Phase = Phases[ii]
    dPhase = Phase[1] - Phase[0]
    if dPhase < 0:
        dPhase += 360
    dPhase = dPhase * 2 * np.pi / 360 # convert to radians
    mks = 3*np.log(Amp[0]/5)+2
    ax.plot(Om/om0, dAmp, '*r', markersize=mks)
    ax2.plot(Om/om0, dPhase/np.pi, '*b', markersize=mks)
    ii += 1
    
plt.show()
