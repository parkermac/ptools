"""
Code to find out how many "winners" there are in a given particulator
experiment, and compare it to a driver like the AB8d (Austin Barth 8 day)
wind stress.
"""

# USER: set values
pmdir = '/Users/PM3/Documents/'
moor = 'C2005_RN'
matfile = pmdir + 'tools_output/moor_out/' + moor + '.mat'

pdir = pmdir + 'tools_output/particulator_out/jdf_inflow_10d/' 

testing = True

# IMPORTS
import matplotlib.pyplot as plt
import numpy as np
import zfun; reload(zfun)
import matfun; reload(matfun)

# Mooring record - load and parse
matdata = matfun.loadmat(matfile)
M = matdata['M']
td = M['td']
mtd = td - td[0] # time in days from the start of the year
varname = 'svstr'
tau = M[varname]
smooth = zfun.filt_AB8d(tau)

# particulator data
import netCDF4 as nc
    
win_list = []

yds = np.arange(5, 350 + 5, 5, dtype=int)
wins = np.nan * yds
cc = 0
for yd in yds:
    pfn = pdir + 'yd.' + str(yd) + '.nc'
    ds = nc.Dataset(pfn)
    p = dict()
    for vn in ds.variables:
        # print ds.variables[vn] # to see info
        p[vn] =  ds.variables[vn][:]
    ds.close()
    
    ptd = p['t'][:,0]/86400. # time in days from start of year
    
    px = p['lon']
    py = p['lat']
    NT, NP = px.shape
    
    # determine distance from a point and then use that to decide
    # winners and losers
    # center point of circle (from the matlab code)
    lon0 = -1.246791925614623e+02;
    lat0 = 48.482319275124404;
    mlr = np.pi*np.mean(py[0,0])/180.;
    clat = np.cos(mlr);
    earth_rad = 6371e3 # m
    xrad = np.pi * px /180
    yrad = np.pi * py / 180
    x = earth_rad * clat * np.pi * (px - lon0) / 180.
    y = earth_rad * np.pi * (py - lat0) / 180.
    dist = np.sqrt(x**2 + y**2)
    
    win_count = 0
    for ii in range(NP):
        if np.nanmin(dist[:, ii]) < 10e3:
            iend = np.where(dist[:, ii] < 10e3)[0][0]
            if p['z'][iend, ii] < -70.0:
                win_count += 1
    wins[cc] = win_count
    
    cc += 1

# plot the results

yds = yds + 5
    
plt.close()

fig = plt.figure(figsize=(12,8))

ax1 = fig.add_subplot(211)

ax1.plot(mtd, smooth, '-k')
ax1.fill_between(mtd, smooth, where=(smooth<0), color='aqua')#, alpha=.1)
ax1.fill_between(mtd, smooth, where=(smooth>0), color='tomato')
ax1.set_xlim([0, 365])
ax1.grid()
ax1.axhline(0)
#ax1.axhspan(0, .1)

ax2 = ax1.twinx()
ax2.plot(yds, wins, 'or')
ax2.set_ylabel('Wins', color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')
ax2.set_xlim([0, 365])
    

#ax.set_xlabel('Time (days)')
ax1.set_ylabel('Low-passed N-S Windstress (Pa)')
ax1.set_title(moor + ' ' + varname)


ax = fig.add_subplot(212)
smi = np.interp(yds, mtd, smooth) # smoothed wind at times yds
wdown = np.mean(wins[np.where(smi>0)])
wup = np.mean(wins[np.where(smi<=0)])
ax.plot(smi[smi>0], wins[smi>0], 'o', color=(.9,.5,.5))
ax.plot(smi[smi<0], wins[smi<0], 'oc')
ax.grid()
ax.set_xlabel('smooth')
ax.set_ylabel('Number of Wins')

ax.text(.1, .9, 'Mean during upwelling ' + str(round(wup, 1)), transform=ax.transAxes)
ax.text(.9, .9, 'Mean during downwelling ' + str(round(wdown, 1)), horizontalalignment='right', transform=ax.transAxes)


plt.show()
