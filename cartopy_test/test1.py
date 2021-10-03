"""
From Rob H.
"""
import numpy as np
import matplotlib.pyplot as plt

import cartopy.feature as cfeature
cl = cfeature.GSHHSFeature('high') # 'full'
geom_list = [item for item in cl.geometries()]
g = geom_list[7] #australia

"""
This threw a warning:

/Applications/miniconda3/envs/loenv/lib/python3.9/site-packages/shapefile.py:391:
UserWarning: Shapefile shape has invalid polygon:
no exterior rings found (must have clockwise orientation);
interpreting holes as exteriors.
"""

# a.exterior.xy is a tuple of points
x = g.exterior.xy[0]
y = g.exterior.xy[1]

plt.close('all')
fig = plt.figure()
ax = fig.add_subplot(111)

#ax.plot(*a.exterior.xy, '-g')
# the * is a shorthand for unpacking the tuple - too cute for me
ax.plot(x,y,'-b')
plt.show()
