"""
Code to test argrelextrema.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

x = np.linspace(-6,6,1000)
y = -x**2 + np.sin(x)

Imax = argrelextrema(y, np.greater_equal, order=10)
imax = Imax[0]

plt.close('all')
fig = plt.figure(figsize=(13,8))
ax = fig.add_subplot(111)

ax.plot(x,y, '-c')

if len(imax) > 0:
    for ii in imax:
        ax.plot(x[ii], y[ii], '*')
        
ax.set_title('len(imax) = %d' % (len(imax)))

plt.show()
