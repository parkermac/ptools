"""
Experiments with cartopy.
"""

import cartopy.crs as ccrs
import matplotlib.pyplot as plt

plt.close('all')

ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines(resolution='50m')

"""
This shows a warning the FIRST time you run it with a given resolution, e.g.:

/Applications/miniconda3/envs/loenv/lib/python3.9/site-packages/cartopy/io/__init__.py:241:
 DownloadWarning: Downloading: https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_coastline.zip
  warnings.warn(f'Downloading: {url}', DownloadWarning)
"""

# Save the plot by calling plt.savefig() BEFORE plt.show()
# plt.savefig('coastlines.pdf')
# plt.savefig('coastlines.png')

plt.show()
