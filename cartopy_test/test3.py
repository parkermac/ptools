"""
Based on an example at:
https://cartopy-pelson.readthedocs.io/en/readthedocs/matplotlib/feature_interface.html
"""

import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature


plt.close('all')

ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-130, -122, 42, 52])

# Put a background image on for nice sea rendering.
ax.stock_img()

# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')

ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(states_provinces, edgecolor='gray')

plt.show()
