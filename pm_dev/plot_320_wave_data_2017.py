# -*- coding: utf-8 -*-
"""
Code to plot the wave data gathered by the Ocean 320 class
on January 11, 2017, Portage Bay, and Harris Lab

"""

import numpy as np
import matplotlib.pyplot as plt

# wavelenght [m]
L = np.array([.4, 1.5, .3, .85, 1, .75, .75, .5, .16,
    4.1, 1.46, 1.46, 1.5, 1, .5, .83, .58, 1])

# period [s]
T = np.array([.4, 1.5, 3, 1.29, .83, .78, .8, .52, .15,
    1.2, 2.48, 1.618, 1.5, 1.49, 1, 1.18, 1.35, 1.23])

# water depth [m] values = 0.75 are in Harris
H = np.array([np.nan, np.nan, .2, .9, .75, .28, .3, .3, .3,
    2.9, .75, .75, .2, .2, .2, np.nan, np.nan, np.nan])

# mask for waves that we expect to follow deep water theory
mask = L <= 2*H

# wave speed [ m s-1]
c = L/T

# theoretical curve
g = 9.8 # gravity [m s-2]
LL = np.linspace(0,4.5,1000)
c_dw = np.sqrt(g*LL/(2*np.pi))
c_sw = np.sqrt(g*H)

c_full = np.sqrt( (g*L/(2*np.pi)) * np.tanh(2*np.pi*H/L) )

plt.close('all')

fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot(1,1,1)
ax.plot(L[mask], c[mask], '*r', markersize=16)
ax.plot(LL, c_dw, '-r', linewidth=3)
ax.plot(L[~mask], c[~mask],'ob', markersize=10)
ax.plot(L[~mask], c_full[~mask], 'oy', markersize=10)
ax.legend(('Observations','Deep Water Theory',
           'Maybe water too shallow', 'Full Theory'), loc=2)
ax.grid()
ax.set_xlim(0, 4.5)
ax.set_ylim(0,3.5)
ax.set_xlabel('Wavelength [m]')
ax.set_ylabel('Phase Speed [m s-1]')
ax.set_title('OCN 320 Portage Bay + Harris Lab Observations 1/11/2017')

plt.show()



