# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 16:33:24 2016

@author: PM5

Numerical integration of the 1D shallow water equations.

"""

import numpy as np
import matplotlib.pyplot as plt

g = 9.8
H = 10
a = 1
c = np.sqrt(g*H)
U = c*a/H
Cd = 3e-3
R = Cd*U/H

dx = 100
dt = 1
NX = 100
NT = 5000

L = NX*dx
T = NT*dt

u = np.nan * np.ones((NT, NX))
e = np.nan * np.ones((NT, NX))

ph = np.linspace(0, np.pi, NX)
x = np.linspace(-L/2, L/2, NX)
t = np.linspace(0, (NT-1)*dt, NT)

e[0, :] = a * np.exp(-(x/1000)**2)
u[0, :] = c * (a/H) * np.exp(-(x/1000)**2)

for ii in range(1,NT):
    if ii == 1:
        # forward for first time step
        u[ii, 1:-1] = u[ii-1, 1:-1] + dt*( -g*(e[ii-1, 2:] - e[ii-1, :-2])/dx
                    - R*u[ii-1, 1:-1])
        e[ii, 1:-1] = e[ii-1, 1:-1] -dt*H*(u[ii-1, 2:] - u[ii-1, :-2])/dx
    else:
        # leapfrog thereafter
        u[ii, 1:-1] = u[ii-2, 1:-1] + 2*dt*( -g*(e[ii-1, 2:] - e[ii-1, :-2])/dx
                    - R*u[ii-1, 1:-1])
        e[ii, 1:-1] = e[ii-2, 1:-1] -2*dt*H*(u[ii-1, 2:] - u[ii-1, :-2])/dx

    u[ii, 0] = 0
    u[ii, -1] = 0
    e[ii, 0] = e[ii, 1]
    e[ii, -1] = e[ii, -2]

# plotting

plt.close()

fig = plt.figure(figsize=(15,10))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

tfac = 100
offset = np.linspace(0,T/tfac,NT).reshape(NT,1) * np.ones(NX).reshape(1,NX)

nlines = 100
step = int(NT/nlines)
uu = u[::step, :].T + offset[::step, :].T
ee = e[::step, :].T + offset[::step, :].T

ax1.plot(x,uu)
ax2.plot(x,ee)

ax1.set_title('u')
ax2.set_title('eta')

ax1.set_xlabel('x (m)')
ax2.set_xlabel('x (m)')
ax1.set_ylabel('t/' + str(tfac) + ' (s)')

ax1.set_xlim(-L/2, L/2)
ax1.set_ylim(0, T/tfac)
ax2.set_xlim(-L/2, L/2)
ax2.set_ylim(0, T/tfac)

plt.show()
