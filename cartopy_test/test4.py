"""
Based on this page from Rob H.:

https://github.com/hetland/python4geosciences/blob/master/materials/7_shapefiles.ipynb

RESULT: The highest resolution gshhs coastline ('full') is identical to the product
I used for my Canadian coastline, and appears to be a bit different and in my opinion
not as nice as the US coastline I use - which I got a long time ago from the coastline
extractor.

Performance for scale='f' and athe cas6 domain:
1517 rings found
Took 5.19 sec
- although it seems to slow down when run repeatedly
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.io.shapereader as shpreader
import shapely.geometry
from time import time
import xarray as xr

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

# load the new Puget Sound bathymetry
dir0 = Ldir['data'] / 'topo' / 'psdem_10m'
fn = dir0 / 'puget_sound_13_navd88_2014.nc'
ds = xr.open_dataset(fn)

step = 5 # downsample
x = ds.lon.values[::step]
y = ds.lat.values[::step]
z = ds.Band1.values[::step,::step]

xx, yy = np.meshgrid(x,y)
xp, yp = pfun.get_plon_plat(xx,yy)


g_shp = shpreader.gshhs(scale='f')
# Scales: ‘auto’, ‘coarse’, ‘low’, ‘intermediate’, ‘high, or ‘full’ (default is ‘auto’),
# or the first letter of the choice.
reader = shpreader.Reader(g_shp)
recs = reader.records()

aa = (-130, -122, 42, 52)
x0, x1, y0, y1 = aa

tt0 = time()
cl_list = []
for cl in recs:
    """
    Each cl is a "Record": "A single logical entry from a shapefile,
    combining the attributes with their associated geometry."
    
    In other words, it is a shapely.geometry.polygon.Polygon object
    and some extra information.
    
    It has: attributes, bounds, and geometry.
    
    cl.attributes is a dict, but these are not too informative for gshhg
    
    cl.bounds is a tuple like (-9.500389, 1.2695, 180.0, 77.719583)
    which I think is lon,lat of the lower left and upper right corners
    with lon in -180:180 convention
    
    You can get the underlying lon, lat vectors like this:
    x = cl.geometry.exterior.xy[0]
    y = cl.geometry.exterior.xy[1]
    
    cl.geometry.exterior is an instance of a
    shapely.geometry.polygon.LinearRing.
    
    cl.geometry is a shapely.geometry.polygon
    
    """
    xx0, yy0, xx1, yy1 = cl.bounds
    
    if (
        # this gets all the islands
        ((xx0>x0)and(xx0<x1))and((yy0>y0)and(yy0<y1))
        or ((xx1>x0)and(xx1<x1))and((yy0>y0)and(yy0<y1))
        or ((xx0>x0)and(xx0<x1))and((yy1>y0)and(yy1<y1))
        or ((xx1>x0)and(xx1<x1))and((yy1>y0)and(yy1<y1))
        # and this gets the US mainland coastline
        or ((x0>xx0)and(x0<xx1))and((y0>yy0)and(y0<yy1))
        or ((x1>xx0)and(x1<xx1))and((y0>yy0)and(y0<yy1))
        or ((x0>xx0)and(x0<xx1))and((y1>yy0)and(y1<yy1))
        or ((x1>xx0)and(x1<xx1))and((y1>yy0)and(y1<yy1))
        ):
        c_ring = cl.geometry.exterior
        cl_list.append(c_ring)
print('%d rings found' % (len(cl_list)))
print('Took %0.2f sec' % (time()-tt0))

# PLOTTING
plt.close('all')
fig = plt.figure(figsize=(10,13))
ax = fig.add_subplot(111)

cs = ax.pcolormesh(xp,yp,z, vmin=-50, vmax=10, cmap='terrain')
ax.contour(xx,yy,z,[0], colors='r')
fig.colorbar(cs, ax=ax)

for C in cl_list:
    x = C.xy[0]
    y = C.xy[1]
    ax.plot(x,y,'-b')
ax.axis(aa)
pfun.add_coast(ax)
pfun.dar(ax)
plt.show()
