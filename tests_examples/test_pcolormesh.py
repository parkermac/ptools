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

fig = plt.figure(figsize=(20,8))

ax = fig.add_subplot(131)
cs = ax.pcolormesh(xx,yy,zz)
fig.colorbar(cs, ax=ax, orientation='horizontal')
# Result: this aligns the LOWER LEFT corner of each colored square
# with the axes, and DROPS the last row and column of the data.

ax = fig.add_subplot(132)
ax.pcolormesh(xx,yy,zz[:-1, :-1])
fig.colorbar(cs, ax=ax, orientation='horizontal')
# Result: Here we have axes with size one larger than the data, and
# the colored squares are plotted as if the axes referred to the CORNERS
# of each square.  The plot is IDENTICAL to that of ax0 - so it is not doing
# bilinear interpolation.

ax = fig.add_subplot(133)
ax.pcolormesh(xx,yy,zz, shading='gouraud')
fig.colorbar(cs, ax=ax, orientation='horizontal')
# Result: I'm not sure what this does.

# A question to answer is how best to plot masked arrays.

plt.show()