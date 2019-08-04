"""
Code to explore shallow water frictional tides in bays, like
the quarter-wave resonance problem.
"""

import numpy as np
import matplotlib.pyplot as plt

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

# set parameters
g = 9.8 # gravity [m/s2]
H = 20 # channel depth [m]
Cd = 3e-3 # drag coefficient []
om = 1.4e-4 # frequency [s-1] typical of M2 tide
Ut = 1 # scale of tidal velocity [m/s]
L = 500e3 # channel length [m]
a = 1 # tidal amplitude at mouth [m]

# derived parameters
R = Cd*Ut/H
c = np.sqrt(g*H)

x = np.linspace(0, L, 500) # x-coordinate [m]

def get_K(om, c, R):
    K = (om/c) * np.sqrt( np.complex(1, R/om) )
    return K

# solution form for a given om (except for time dependence)
K = get_K(om, c, R)
E = (a * np.exp(1j*K*x))

# PLOTTING
plt.close('all')
fig = plt.figure(figsize=(11,7))

ax = fig.add_subplot(111)
for omt in np.linspace(0, 2*np.pi, 12):
    ax.plot(x/1e3, np.real(E*np.exp(-1j*omt)))
ax.set_xlim(0, L/1e3)
ax.set_xlabel('Along Channel Distance from Mouth (km)')
ax.set_ylabel('Surface Height (m)')


# print some scale information
ax.text(.95, .9, 'H = %d m, R/omega = %0.1f' % (H, R/om),
    transform=ax.transAxes, horizontalalignment='right')
ax.text(.95, .8, 'Tidal Lengthscale = 1/k_real = %d km' % (1/(K.real*1e3)),
    transform=ax.transAxes, horizontalalignment='right')
ax.text(.95, .7, 'Frictional e-Folding Scale = %d km' % (1/(K.imag*1e3)),
    transform=ax.transAxes, horizontalalignment='right')
ax.text(.05, .9, 'Surface height at different times for frequency = M2',
    transform=ax.transAxes, horizontalalignment='left',color='r',fontweight='bold')
ax.grid(True)
    


plt.show()

