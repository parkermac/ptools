"""
Following example from Rob:
https://github.com/hetland/python4geosciences/blob/master/materials/ST_shapes.ipynb
"""
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline

import shapely.geometry
import shapely.ops
import cartopy
import cartopy.io.shapereader as shpreader

point = shapely.geometry.Point(0.2, 1.0)

plt.close('all')
x,y = point.xy
plt.plot(x[0],y[0],'*r')
plt.show()