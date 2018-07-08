"""
Test of nearest-neighbor extrapolation to fill missing (masked) values
on a plaid grid.
"""

import numpy as np
import numpy.ma as ma

# create a data array, with masked areas
a = np.arange(100).reshape(10,10)
fill_value=-99
a[2:4,3:8] = fill_value
a[7:,7:] = fill_value
a = ma.masked_array(a,a==fill_value)

# create axes
xx = np.linspace(0, 100, a.shape[1])
yy = np.linspace(0, 10, a.shape[0])
x, y = np.meshgrid(xx, yy)

# do the extrapolation
from scipy.spatial import cKDTree
xygood = np.array((x[~a.mask],y[~a.mask])).T
xybad = np.array((x[a.mask],y[a.mask])).T
b = a.copy()
b[a.mask] = a[~a.mask][cKDTree(xygood).query(xybad)[1]]
# asking for [1] gives the index of the nearest good value
# and [0] would return the distance, I think.

print(a)

print(b)