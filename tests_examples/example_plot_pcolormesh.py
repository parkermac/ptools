"""
Test of cell alignment in pcolor-like plotting.
"""

import numpy as np
import matplotlib.pyplot as plt

x = np.arange(4) # 0 to 3
y = np.arange(3) # 0 to 2
xx, yy = np.meshgrid(x, y)
zz = xx**2 + yy**2

# plotting
plt.close()

fig = plt.figure(figsize=(20,5))

ax0 = fig.add_subplot(131)
ax0.pcolormesh(xx,yy,zz)
# Result: this aligns the LOWER LEFT corner of each colored square
# with the axes, and DROPS the last row and column of the data.

ax1 = fig.add_subplot(132)
ax1.pcolormesh(xx,yy,zz[:-1, :-1])
# Result: Here we have axes with size one larger than the data, and
# the colored squares are plotted as if the axes referred to the CORNERS
# of each square.  The plot is IDENTICAL to that of ax0 - sot it is not doing
# bilinear interpolation.

ax2 = fig.add_subplot(133)
ax2.pcolormesh(xx,yy,zz, shading='gouraud')
# Result: I'm not sure what this does.

# A question to answer is how best to plot masked arrays.

plt.show()