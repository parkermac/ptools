"""
Code to test control of colorbars
"""

import matplotlib.pyplot as plt
import numpy as np

plt.close('all')
fig = plt.figure()

x = np.linspace(-10,10,100)
y = np.linspace(-10,10,100)
X, Y = np.meshgrid(x,y)
Z = X**2 - Y**2

fs = 16

ax = fig.add_subplot(111)

cs = ax.pcolormesh(X, Y, Z, cmap='Spectral')
ax.tick_params(labelsize=fs)

cb = fig.colorbar(cs)
cb.ax.tick_params(labelsize=fs)

plt.show()
