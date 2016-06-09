# -*- coding: utf-8 -*-
"""
Code to plot the wave data gathered by the Ocean 320 class
on April6, 2016, Portag Bay
Created on Wed Apr  6 16:08:16 2016

@author: PM5
"""

import numpy as np
import matplotlib.pyplot as plt

# wavelenght [m]
L = np.array([1, .5, 2, .3, .25, 2, 3])

# period [s]
T = np.array([.6, 1, 1.5, 2, .5, 1, 1.5])

# water depth [m]
H = np.array([.4, 1, 0.5, .19, 1, 3.3, 1.9])

# mask for waves that we expect to follow deep water theory
mask = L/2 <= H

# wave speed [ m s-1]
c = L/T

# theoretical curve
g = 9.8 # gravity [m s-2]
LL = np.linspace(0,4,1000)
c_dw = np.sqrt(g*LL/(2*np.pi))
c_sw = np.sqrt(g*H)

plt.close('all')

fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot(1,1,1)
ax.plot(L[mask], c[mask], '*r', markersize=16)
ax.plot(LL, c_dw, '-r', linewidth=3)
ax.plot(L[~mask], c[~mask],'ob', markersize=10)
ax.plot(L[~mask], c_sw[~mask],'oc', markersize=10)
ax.legend(('Observations','Deep Water Theory',
           'Maybe water too shallow', 'Shallow Water Theory'))
ax.grid()
ax.set_xlim(0, 4)
ax.set_ylim(0,3)
ax.set_xlabel('Wavelength [m]')
ax.set_ylabel('Phase Speed [m s-1]')
ax.set_title('OCN 320 Portage Bay Observations 4/6/2016')

plt.show()



