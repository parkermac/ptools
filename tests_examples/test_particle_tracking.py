# -*- coding: utf-8 -*-
"""
Code to experiment with particle tracking numerical schemes.

RESULT:

All results are assuming dt = 3600.

As vort increases above 1e-4 the errors of the rk2 methods become evident.

All methods do quite well with just the inertial circle

rk2 and rk2t are very similar.

rk4 and rk4t are very similar.

Never use euler.

Created on Mon Jan 25 12:52:09 2016

@author: PM5
"""

import numpy as np
import matplotlib.pyplot as plt

#%% Flow field

# Assume that t is time in seconds, and x,y are m

# solid body rotation with vorticity "vort"
# which will complete a circle in time 2*pi/(vort/2)
# or 8.7 hours for vort = 1e-4
#
# plus an inertial oscillation (period 17.5h for f=1e-4) with velocity "U"
# which on its own will cause particles to go in a cyclonic circle of
# diameter U*T/pi where T = 17.5 hours (4 km for U = 0.2)

vort = 0*2e-4 # rad s-1
U = .4 # m s-1
 
def u(x,y,t):
    return -(vort/2)*y + U*np.cos(2*np.pi*t/(3600*18))
def v(x,y,t):
    return (vort/2)*x + U*np.sin(2*np.pi*t/(3600*18))

dt = 3600
ndays = 10   
t = np.arange(0,86400*ndays,dt)
n = len(t)
x = np.nan * t
y = np.nan * t

x[0] = 1000
y[0] = 0

#%% define integration methods

def track(x, y, t, dt, method):
    
    if method == 'euler':
        ii = 0
        for tt in t[1:]:
            ii += 1       
            # Forward Euler
            x0 = x[ii-1]
            y0 = y[ii-1]
            t0 = t[ii-1]            
            u0 = u(x0, y0, t0)
            v0 = v(x0, y0, t0)
            
            x[ii] = x0 + dt*u0
            y[ii] = y0 + dt*v0
               
    elif method == 'rk2':
        ii = 0
        for tt in t[1:]:
            ii += 1
            # RK2, or midpoint method, or leapfrog
            x0 = x[ii-1]
            y0 = y[ii-1]
            t0 = t[ii-1]            
            u0 = u(x0, y0, t0)
            v0 = v(x0, y0, t0)
            
            x1 = x0 + dt*u0/2
            y1 = y0 + dt*v0/2
            thalf = t0 + dt/2
            u1 = u(x1, y1, thalf)
            v1 = v(x1, y1, thalf) 
            
            x[ii] = x0 + dt*u1
            y[ii] = y0 + dt*v1
            
    elif method == 'rk2t':
        ii = 0
        for tt in t[1:]:
            ii += 1
            # RK2, or midpoint method, or leapfrog
            # but assuming we only have velocity at discrete times (dt)
            x0 = x[ii-1]
            y0 = y[ii-1]
            t0 = t[ii-1]            
            u0 = u(x0, y0, t0)
            v0 = v(x0, y0, t0)
            
            x1 = x0 + dt*u0/2
            y1 = y0 + dt*v0/2
            tfull = t0 + dt
            u1 = (u(x1, y1, t0) + u(x1, y1, tfull))/2
            v1 = (v(x1, y1, t0) + v(x1, y1, tfull))/2 
            
            x[ii] = x0 + dt*u1
            y[ii] = y0 + dt*v1 
              
    elif method == 'rk4':
        ii = 0
        for tt in t[1:]:
            ii += 1
            # RK4
            x0 = x[ii-1]
            y0 = y[ii-1]
            t0 = t[ii-1]            
            u0 = u(x0, y0, t0)
            v0 = v(x0, y0, t0)
            
            x1 = x0 + dt*u0/2
            y1 = y0 + dt*v0/2
            thalf = t0 + dt/2
            u1 = u(x1, y1, thalf)
            v1 = v(x1, y1, thalf) 
            
            x2 = x0 + dt*u1/2
            y2 = y0 + dt*v1/2
            u2 = u(x2, y2, thalf)
            v2 = v(x2, y2, thalf) 
            
            x3 = x0 + dt*u2
            y3 = y0 + dt*v2
            tfull = t0 + dt            
            u3 = u(x3, y3, tfull)
            v3 = v(x3, y3, tfull)
            
            x[ii] = x0 + dt*(u0 + 2*u1 + 2*u2 + u3)/6
            y[ii] = y0 + dt*(v0 + 2*v1 + 2*v2 + v3)/6

    elif method == 'rk4t':
        ii = 0
        for tt in t[1:]:
            ii += 1
            # RK4
            # but assuming we only have velocity at discrete times (dt)
            # requires 6 velocity extractions, versus 4 for rk4,
            # and 3 for rk2t and 2 for rk2, and 1 for euler
            x0 = x[ii-1]
            y0 = y[ii-1]
            t0 = t[ii-1]            
            u0 = u(x0, y0, t0)
            v0 = v(x0, y0, t0)
            
            x1 = x0 + dt*u0/2
            y1 = y0 + dt*v0/2
            tfull = t0 + dt
            u1 = (u(x1, y1, t0) + u(x1, y1, tfull))/2
            v1 = (v(x1, y1, t0) + v(x1, y1, tfull))/2 
            
            x2 = x0 + dt*u1/2
            y2 = y0 + dt*v1/2
            u2 = (u(x2, y2, t0) + u(x2, y2, tfull))/2
            v2 = (v(x2, y2, t0) + v(x2, y2, tfull))/2 
            
            x3 = x0 + dt*u2
            y3 = y0 + dt*v2
            tfull = t0 + dt            
            u3 = u(x3, y3, tfull)
            v3 = v(x3, y3, tfull)
            
            x[ii] = x0 + dt*(u0 + 2*u1 + 2*u2 + u3)/6
            y[ii] = y0 + dt*(v0 + 2*v1 + 2*v2 + v3)/6
            
    return x, y

#%% plotting   
plt.close()
fig = plt.figure()

ax = fig.add_subplot(111)

clist = ['r','b','g','c','m','k']
ytlist = [.95,.9,.85,.8,.75,.7]

counter = 0
for method in ['rk2', 'rk2t', 'rk4', 'rk4t']:
    x, y = track(x, y, t, dt, method)
    ax.plot(x/1000, y/1000, '-', color=clist[counter])
    ax.text(.05, ytlist[counter], method, fontsize=16,
            color=clist[counter],
            transform=ax.transAxes)
    counter += 1
        
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')

#ax.set_title(method)
ax.grid()
ax.axis('square')

plt.show()
