# importing of NetCDF files
import sys, os


import numpy as np

import netCDF4 as nc


sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
sys.path.append(os.path.abspath('../../LiveOcean/plotting'))
import pfun
ds = nc.Dataset('../../LiveOcean_data/grids/cas6/grid.nc')
latbounds = [47.581, 47.663]
lonbounds = [-122.441, -122.331]

# PM: use lon_psi and lat_psi as the corners for pcolormesh
latp = ds['lat_psi'][:].data
lonp = ds['lon_psi'][:].data
# PM: but use lon_rho and lat_rho (grid cell centers) for your release points
lats = ds['lat_rho'][:].data
lons = ds['lon_rho'][:].data
maskr = ds['mask_rho'][:].data
# PM: get h field for plotting
h = ds['h'][:].data
h[maskr==0] = np.nan

        # latitude lower and upper index

# PM: I changed lats[:,1] to lats[:,0], and etc. because python
# starts at zero (makes no differnce because the arrays are plaid)
rlat = lats[:,0]>=latbounds[0]
latin = np.argmin(~rlat)

rlat = lats[:,0]<=latbounds[1]
latax = np.argmax(~rlat)

rlon = lons[0,:]>=lonbounds[0]
lonin = np.argmin(~rlon)

rlon = lons[0,:]<=lonbounds[1]
lonax = np.argmax(~rlon)

maskr = maskr[latin:latax, lonin:lonax]
print(maskr)
lat = lats[latin:latax, 0]
print(len(lat))
lon = lons[0, lonin:lonax]

      # 1=water, 0=land
il=[]
jl=[]
paper=np.zeros_like(maskr)
for i in range(1,15):
    for j in range(1,15):
        if maskr[i, j]==1 and maskr[i-1, j]==0:
                paper[i,j]=2
                il.append(i)
                jl.append(j)
        elif maskr[i, j]==1 and maskr[i+1, j]==0:
                paper[i,j]=2
                il.append(i)
                jl.append(j)
        elif maskr[i, j]==1 and maskr[i, j+1]==0:
                paper[i,j]=2
                il.append(i)
                jl.append(j)
        elif maskr[i, j]==1 and maskr[i, j-1]==0:
                paper[i,j]=2
                il.append(i)
                jl.append(j)
        elif maskr[i, j]==1 and maskr[i-1, j-1]==0:
                paper[i,j]=2
                il.append(i)
                jl.append(j)
        elif maskr[i, j]==1 and maskr[i-1, j+1]==0:
                paper[i,j]=2
                il.append(i)
                jl.append(j)
        elif maskr[i, j]==1 and maskr[i+1, j-1]==0:
                paper[i,j]=2
                il.append(i)
                jl.append(j)
        elif maskr[i, j]==1 and maskr[i+1, j+1]==0:
                paper[i,j]=2
                il.append(i)
                jl.append(j)
        else:
            continue
print(paper)
latf = lat[il]
lonf = lon[jl]

# PM: PLOTTING to check your results
import matplotlib.pyplot as plt
plt.close('all')
fig = plt.figure(figsize=(10,12))
ax = fig.add_subplot(111)

ax.pcolormesh(lonp, latp, h[1:-1, 1:-1], cmap='jet', vmin=0, vmax=200)

pfun.draw_box(ax, lonbounds + latbounds, linewidth=3, color='g')

pfun.dar(ax)
pfun.add_coast(ax)
ax.axis([-122.5, -122.3, 47.5, 47.7])

ax.plot(lonf, latf, '.k')

plt.show()
