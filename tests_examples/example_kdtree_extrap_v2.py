"""
Test of interpolation using cKDTree

Rob says there are problems with this, and he prefers Delaunay.
"""

from scipy.spatial import cKDTree

import numpy as np

x = np.linspace(0, 10, 6)
y = np.linspace(0, 10, 6)
X, Y = np.meshgrid(x, y)

Z = X**2 + Y**2

# note that the Transpose is needed because of what
# cKDTree expects
XY = np.array((X.flatten(), Y.flatten())).T
XYI = np.array((2., 2. )).T

ZI = Z.flatten()[cKDTree(XY).query(XYI)[1]]

# Here is how you would use this for getting nearest neighbors
# on a 1D axis.  Just use arrays of single-element tulples
a = np.array((x.flatten(), )).T
ai = np.array((3., )).T
# ai = (3.,) # also works for the trivial case of one point
# b[0] is the distance to the nearest neighbor, and b[1] is the index
b = cKDTree(a).query(ai)