"""
Plot results of LO_hycom_to_netcdf.py.
"""
import zfun; reload(zfun)
nc_dir = '../../tools_output/pydev_out/hycom_test/'

import netCDF4 as nc

# for reference, here are the choices for fld_name
# fld_list = ['ssh', 's3d', 't3d', 'u3d', 'v3d']
fld_name = 'v3d'

fn = nc_dir + fld_name + '.nc'
ds = nc.Dataset(fn)
fld = ds.variables[fld_name][:]
lon = ds.variables['lon'][:]
lat = ds.variables['lat'][:]
dt = ds.variables['dt'][:]

# make a smoothed version of fld, to remove intertial oscillations
flds = fld.copy()
if fld_name == 'ssh':    
    flds[2:-2,:,:] = (fld[:-4,:,:] + fld[1:-3,:,:] + fld[2:-2,:,:] + fld[3:-1,:,:] + fld[4:,:,:])/5
    flds[:2,:,:] = flds[2,:,:]
    flds[-2:,:,:] = flds[-3,:,:]
else:
    flds[2:-2,:,:,:] = (fld[:-4,:,:,:] + fld[1:-3,:,:,:] + fld[2:-2,:,:,:] + fld[3:-1,:,:,:] + fld[4:,:,:,:])/5
    flds[:2,:,:,:] = flds[2,:,:,:]
    flds[-2:,:,:,:] = flds[-3,:,:,:]

# get time series
nr = 80; nc = 40
if fld_name == 'ssh':
    fld_top = fld[0,:,:].squeeze()    
    fs = fld[:,nr,nc].squeeze()
    fss = flds[:,nr,nc].squeeze()
else:
    nz = -1
    fld_top = fld[0,nz,:,:].squeeze()
    fs = fld[:,nz,nr,nc].squeeze()
    fss = flds[:,nz,nr,nc].squeeze()
    
# get the coastline
coast_file = '../../tools_data/pydev_data/pnw_coast_combined.mat'
import scipy.io
cmat = scipy.io.loadmat(coast_file)
clon = cmat['lon'].squeeze()
clat = cmat['lat'].squeeze()

# PLOTTING
import matplotlib.pyplot as plt
plt.close()

fig = plt.figure(figsize=(15, 8))
cmap = plt.get_cmap(name='rainbow')

ax0 = fig.add_subplot(121)
cs0 = ax0.pcolormesh(lon, lat, fld_top,  cmap = cmap)
ax0.plot(lon[nr,nc], lat[nr,nc],'*r')
# add coastline
ax0.plot(clon, clat, '-k', linewidth=.5)
# limits, scaling, and labels
ax0.axis([lon.min(), lon.max(), lat.min(), lat.max()])
zfun.dar(ax0)
ax0.set_xlabel('Longitude')
ax0.set_ylabel('Latitude')
ax0.set_title(fld_name + ' ' + str(dt[0]/86400.))
fig.colorbar(cs0)

ax1 = fig.add_subplot(222)
ax1.plot((dt - dt[0])/86400, fs,'-*b', linewidth=3 )
ax1.plot((dt - dt[0])/86400, fss,'-g', linewidth=3 )
plt.show()
