"""
This code works out some things related to the Masse (1990) solution for estuarine
inflow.  It is in the limit of a narrow (Dirac delta function) estuary mouth.
"""

import numpy as np
import matplotlib.pyplot as plt


# parameters
f = 1e-4 # Coriolis [s-1]
Q = -1000 # transport [m3/s] negative = inflow
g = 10 # gravity [m s-2]
s = 1e-3 # slope
r = 3e-3 * .2 # Cd * U [m/s]

K = r / (s*f)

y = np.linspace(0,5e4, 1000) # y coordinate [m]

h = s*y # bottom depth [m]

def get_eta(x):
    eta = ((2*f*Q/(g*s))/(np.sqrt(4*np.pi*K*x))) * np.exp((-y**2)/(4*K*x))
    return eta

x = 1e4 # x position to look at [m]

# get the solution for eta at the chosen x
eta = get_eta(x)

# derivatives
dy = y[1] - y[0]
yy = y[:-1] + dy/2 # y (on centers, meaning in-between the standard y points)
etay = np.diff(eta)/dy # deta/dy (on y, meaning on the standard y points)

# make u velocity
UU = (-g/f) * etay # geostrophic u velocity (on centers)
hh = h[:-1] + np.diff(h)/2 # h (on centers)
q = UU*hh # along-coast transport (on centers)
# interpolate to make u on the regular y coordinate
U = np.nan * eta
U[1:-1] = UU[:-1] + np.diff(UU)/2 # (on y)

# make v-velocity (includes ageostrophic part)
etax = (get_eta(x+1e3) - get_eta(x-1e3))/2e3 # deta/dx
V = (g/f)*etax + r*U/(f*h)# (on y)

# make the friction term from the vorticity equation [s-2]
friction = np.nan * eta
Uonh = UU/hh
friction[1:-1] = r * (Uonh[1:] - Uonh[:-1])/dy # (on y)

# make the stretching term from the vorticity equation [s-2]
stretching = f * V * s / h # (on y)

# make relative vorticity (approximate as -du/dy)
etayy = np.nan * eta
etayy[1:-1] = np.diff(etay)/dy # (on y)
Zeta = g * etayy / f # relative vorticity (on y)

# PLOTTING
plt.close('all')

fs=16
plt.rc('font',size=fs)

fig = plt.figure(figsize=(12,12))

NR = 5

ax = fig.add_subplot(NR,1,1)
ax.plot(y/1e3,eta*100)
ax.set_ylabel('eta [cm]')
ax.set_title('%d km downstream of mouth' % (x/1e3))
ax.grid(True)
ax.set_xticklabels([])
ax.set_xlim(0, y.max()/1e3)

ax = fig.add_subplot(NR,1,2)
ax.plot(yy/1e3, q)
ax.set_ylabel('h * U')
ax.grid(True)
ax.text(.95, .05, 'Q = %d [m3/s]' % (q.sum()*dy), ha='right', transform=ax.transAxes, weight='bold')
ax.set_xticklabels([])
ax.set_xlim(0, y.max()/1e3)

ax = fig.add_subplot(NR,1,3)
ax.plot(y/1e3, Zeta/f)
ax.set_ylabel('Zeta/f')
ax.axhline()
ax.grid(True)
ax.set_xticklabels([])
ax.set_xlim(0, y.max()/1e3)

ax = fig.add_subplot(NR,1,4)
# test of Masse balance.  RESULT: looks good
ax.plot(y/1e3, etax, '--g', lw=3)
ax.plot(y/1e3, K*etayy, '-m', lw=2) 
ax.text(.95, .85, 'deta/dx', c='g', ha='right', transform=ax.transAxes, weight='bold')
ax.text(.95, .7, 'K*d2eta/dy2', c='m', ha='right',  transform=ax.transAxes, weight='bold')
ax.axhline()
ax.set_xticklabels([])
ax.grid(True)
ax.set_xlim(0, y.max()/1e3)

ax = fig.add_subplot(NR,1,5)
# vorticity equation terms [s-2]
ax.plot(y/1e3, stretching, '-r', y/1e3, friction, '-b')
ax.text(.95, .85, 'Stretching', c='r', ha='right', transform=ax.transAxes, weight='bold')
ax.text(.95, .7, 'Friction', c='b', ha='right',  transform=ax.transAxes, weight='bold')
ax.axhline()
ax.set_xlabel('Y [km]')
ax.grid(True)
ax.set_xlim(0, y.max()/1e3)

plt.show()
plt.rcdefaults()