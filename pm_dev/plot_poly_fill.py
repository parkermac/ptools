"""
Test of filling polygons in matplotlib

"""

import numpy as np
import matplotlib.pyplot as plt

x = np.arange(100)
y = np.arange(100)
X, Y = np.meshgrid(x, y)
Z = np.sin(X*Y/1000)

xy = np.array([[10,10], [20,10], [20,20], [10,20]])

# plotting

plt.close()

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

ax.pcolormesh(X, Y, Z)

if False:
# this works
    from matplotlib.patches import Polygon
    p = Polygon(xy)
    ax.add_patch(p)
else:
    # but this seems easier
    ax.fill(xy[:,0], xy[:,1], facecolor='r', alpha=.3, edgecolor='none')


plt.show()