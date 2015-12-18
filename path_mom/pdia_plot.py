"""
Code to explore the diagnostic balances along tracks.
"""

# USER: set values
pmdir = '/Users/PM3/Documents/'


# IMPORTS
import matplotlib.pyplot as plt
import numpy as np
import zfun; reload(zfun)

# make a list of roms files
rpath = pmdir + 'roms/output/D2005_dia/'
fn = rpath + 'ocean_dia_1800.nc'

# get the roms grid info
G, S = zfun.get_basic_info(fn, getT=False)
lon = G['lon_rho']
lat = G['lat_rho']
zbot = -G['h']

import cPickle
pin = open(pmdir + 'tools_output/pydev_out/p75.pkl', 'rb')
p = cPickle.load(pin)
pin.close()

px = p['lon']
py = p['lat']
NT, NP = px.shape

# determine distance from a point and then use that to decide
# winners and losers
# center point of circle (from the matlab code)
lon0 = -1.246791925614623e+02;
lat0 = 48.482319275124404;
clat = np.cos(np.deg2rad(lat0));
earth_rad = 6371e3 # m
x = earth_rad * clat * np.deg2rad(px - lon0)
y = earth_rad * np.deg2rad(py - lat0)
dist = np.sqrt(x**2 + y**2)

# also work with initial positions, and get angle from lon0, lat0
x0 = x[0, :]
y0 = y[0, :]
z0 = p['z'][0, :]
zb0 = - p['H'][0, :]

r0 = np.sqrt(x0[0]**2 + y0[0]**2) # radius in meters
a0 = np.arctan2(y0, x0) # radians
# make it so that a0 is radians around the circle with zero at pi
# moving CCW is still a positive angular change
mNeg = a0 < 0.
mPos = a0 >= 0.
a0[mNeg] = np.pi + a0[mNeg]
a0[mPos] = a0[mPos] - np.pi
if False:
    deepMask = zb0 == zb0.min()
    idm = np.where(deepMask)[0][0]  
    a00 = a0 - a0[idm] # Optional reset to have zero at deep channel
else:
    a00 = a0
# make a distance vector for circumference around the deployment arc
d00 = a00 * r0 # meters
d00k = d00/1000 # km

# get the coastline
coast_file = pmdir + 'tools_data/geo_data/coast/pnw_coast_combined.mat'
import scipy.io
cmat = scipy.io.loadmat(coast_file)
cstlon = cmat['lon'].squeeze()
cstlat = cmat['lat'].squeeze()

# determine winners
iwin = []
iend_list = []
for ii in range(NP):
    if np.nanmin(dist[:, ii]) < 10e3:
        iend = np.where(dist[:, ii] < 10e3)[0][0]
        iend_list.append(iend) 
        if p['z'][iend, ii] < -70.0:
            iwin.append(ii)                
    else:
        iend_list.append(NT) 

t = p['t'][:, 0]
td = t/86400 # time in days from the start of the year

# diagnostic balances along a path
ts_list = ['rate','xadv','yadv','vadv','hdiff','vdiff']
uv_list = ['accel','xadv','yadv','vadv','cor','prsgrd','vvisc']

pdia = dict()
for ts in ts_list:
    pdia['salt_' + ts] = p['salt_' + ts]
    pdia['temp_' + ts] = p['temp_' + ts]
for uv in uv_list:
    pdia['u_' + uv] = p['u_' + uv]
    pdia['v_' + uv] = p['v_' + uv]

# combinations for informative balances - all as if on RHS

TS = dict()
TS['St'] = - pdia['salt_rate'] + pdia['salt_xadv'] + pdia['salt_yadv'] + pdia['salt_vadv']
TS['Sdiff'] = pdia['salt_hdiff'] + pdia['salt_vdiff'] 

TS['Tt'] = - pdia['temp_rate'] + pdia['temp_xadv'] + pdia['temp_yadv'] + pdia['temp_vadv']
TS['Tdiff'] = pdia['temp_hdiff'] + pdia['temp_vdiff'] 

A = dict()
A['Ut'] = - pdia['u_accel'] + pdia['u_xadv'] + pdia['u_yadv'] + pdia['u_vadv']
A['Ucor'] = pdia['u_cor']
A['Upg'] = pdia['u_prsgrd']
A['Ucpg'] = pdia['u_cor'] + pdia['u_prsgrd']
A['Uvisc'] = pdia['u_vvisc']

A['Vt'] = - pdia['v_accel'] + pdia['v_xadv'] + pdia['v_yadv'] + pdia['v_vadv']
A['Vcor'] = pdia['v_cor']
A['Vpg'] = pdia['v_prsgrd']
A['Vcpg'] = pdia['v_cor'] + pdia['v_prsgrd']
A['Vvisc'] = pdia['v_vvisc']

# filter and average

# still need to cut off length of winners

filtlen = 3*40 # 40 hour Hanning window

TSf = TS.copy()
Af = A.copy()
TSfm = dict()
Afm = dict()
for vname in TS.keys():
    vv_temp = np.nan * np.ones(NP)
    for ii in range(NP):
        vv = TS[vname][:, ii]
        vv[iend_list[ii]:] = np.nan
        TSf[vname][:, ii] = zfun.filt_hanning(vv, filtlen)
        vv_temp[ii] = np.nanmean(vv)
    TSfm[vname] = vv_temp
for vname in A.keys():
    vv_temp = np.nan * np.ones(NP)
    for ii in range(NP):
        vv = A[vname][:, ii]
        vv[iend_list[ii]:] = np.nan
        Af[vname][:, ii] = zfun.filt_hanning(vv, filtlen)
        vv_temp[ii] = np.nanmean(vv)
    Afm[vname] = vv_temp 

# plotting

plt.close()

fig = plt.figure(figsize=(16,8))

if False:
    tag_list = ['t', 'cor', 'pg', 'visc']
    c_list = ['r', 'b', 'k', 'm']
    qscl = 0.01
else:    
    tag_list = ['t', 'cpg', 'visc']
    c_list = ['r', 'c', 'm']
    qscl = 0.001
c_dict = dict(zip(tag_list, c_list))

# winners
ax = fig.add_subplot(121)
plt.rcParams['contour.negative_linestyle'] = 'solid'
ax.contour(lon, lat, zbot, [-2000, -1000, -500, -200, -100], colors='g')
ax.axis([-126.5, -123.5, 47, 49.5])
ax.plot(lon0, lat0,'oy', markersize=10)
# add dynamical arrows
for tag in tag_list:
    for ii in range(NP):
        if ii in iwin:
            U = Afm['U' + tag][ii]
            V = Afm['V' + tag][ii]
            # by default quiver has arrows pointing the right way!
            # even when we call ax.set_aspect (from zfun.dar) to fix the underlying map
            ax.quiver(px[0, ii], py[0, ii], U, V, color=c_dict[tag], units='y', scale=qscl, width=.005)
            ax.plot(px[:, ii], py[:,ii], '-k', linewidth=.2)    
zfun.dar(ax)
# add coastline
ax.plot(cstlon, cstlat, '-k', linewidth=.5)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Winners')

# losers
ax = fig.add_subplot(122)
plt.rcParams['contour.negative_linestyle'] = 'solid'
ax.contour(lon, lat, zbot, [-2000, -1000, -500, -200, -100], colors='g')
ax.axis([-126.5, -123.5, 47, 49.5])
ax.plot(lon0, lat0,'oy', markersize=10)
# add dynamical arrows
for tag in tag_list:
    for ii in range(NP):
        if ii not in iwin:
            U = Afm['U' + tag][ii]
            V = Afm['V' + tag][ii]
            ax.quiver(px[0, ii], py[0, ii], U, V, color=c_dict[tag], units='y', scale=qscl, width=.005)
            ax.plot(px[:, ii], py[:,ii], '-k', alpha=.2, linewidth=.2)    
zfun.dar(ax)
# add coastline
ax.plot(cstlon, cstlat, '-k', linewidth=.5)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Losers')

fig2 = plt.figure(figsize=(15, 15))
ax = fig2.add_subplot(111)

cs0 = p['cs'][0, :] # all initial cs values
css = set(np.round(cs0*100)) # get just unique values
# need to round because values are not exact
nloc = len(css) # number of unique initial cs values
ax.plot(d00k[:NP/nloc], zb0[:NP/nloc], '-', color='0.5', linewidth=3)
ax.plot(d00k, z0, 'oc', markersize=3)
ax.plot(d00k[iwin], z0[iwin], '*r', markersize=15)

qscl2 = qscl/50
for tag in tag_list:
    for ii in range(NP):
        U = Afm['U' + tag][ii]
        V = Afm['V' + tag][ii]
        # rotate to point toward the mouth of JdF
        cal = np.cos(a00[ii])
        sal = np.sin(a00[ii])
        UU = U*cal + V*sal
        VV = V*cal - U*sal
        # plot so pointing up means pointing at JdF
        ax.quiver(d00k[ii], z0[ii], -VV, UU, color=c_dict[tag], units='x', scale=qscl2, width=.5)

plt.show()

