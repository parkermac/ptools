"""
Test of plotting bathymetry with shading to make it look cool.
"""

from pathlib import Path
# path to topo data
pth = Path(__file__).absolute().parent.parent.parent / 'ptools_data' / 'topo'

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource

#ds = nc.Dataset(pth / 'cascadia' / 'cascadia_gridded.nc')
ds = nc.Dataset(pth / 'psdem' / 'psdem_2005_plaid_91m.nc')
zz = ds['elev'][:]
z = zz.data
z[zz.mask] = 0
z[np.isnan(z)] = 0
z = np.flipud(z)

for vn in ds.variables:
    print('%s: %s' % (vn, ds[vn].shape))
    
ls = LightSource(azdeg=315, altdeg=45)
cmap = plt.cm.gist_earth

plt.close('all')
fig = plt.figure(figsize=(10,10))

ax = fig.add_subplot(111)

#ax.imshow(ls.hillshade(z), cmap='gray')

rgb = ls.shade(z, cmap=cmap, blend_mode='overlay')
ax.imshow(rgb)


plt.show()