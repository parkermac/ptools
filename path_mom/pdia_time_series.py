"""
Code to explore the diagnostic balances along tracks.

Plots dynamical balances along a single track.
"""

# USER: set values
pmdir = '/Users/PM3/Documents/'

# IMPORTS
import matplotlib.pyplot as plt
import numpy as np
import zfun; reload(zfun)

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
# make a distance vector for circumference around the deployment arc
# zero is due west of the mouth of JdF
# positive is north of this (moving CW around circle)
d0 = a0 * r0 / 1000 # km

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

# PROCESS DYNAMICAL BALANCES
doTS = False
# diagnostic balances along a path
uv_list = ['accel','xadv','yadv','vadv','cor','prsgrd','vvisc']
if doTS:
    ts_list = ['rate','xadv','yadv','vadv','hdiff','vdiff']

pdia = dict()

for uv in uv_list:
    pdia['u_' + uv] = p['u_' + uv]
    pdia['v_' + uv] = p['v_' + uv]
if doTS:
    for ts in ts_list:
        pdia['salt_' + ts] = p['salt_' + ts]
        pdia['temp_' + ts] = p['temp_' + ts]

# combinations for informative balances - all as if on RHS
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
if doTS:
    TS = dict()
    TS['St'] = - pdia['salt_rate'] + pdia['salt_xadv'] + pdia['salt_yadv'] + pdia['salt_vadv']
    TS['Sdiff'] = pdia['salt_hdiff'] + pdia['salt_vdiff'] 
    TS['Tt'] = - pdia['temp_rate'] + pdia['temp_xadv'] + pdia['temp_yadv'] + pdia['temp_vadv']
    TS['Tdiff'] = pdia['temp_hdiff'] + pdia['temp_vdiff'] 

# filter and average
filtlen = 3*40 # 40 hour Hanning window
# have to use some dilligence in copying A to initialize Af
# I think because np.copy just works on ndarrays
Af = dict()
for vname in A.keys():
    Af[vname] = A[vname].copy()
Afm = dict()
for vname in A.keys():
    vv_temp = np.nan * np.ones(NP)
    for ii in range(NP):
        vv = A[vname][:, ii].copy()
        vv[iend_list[ii]:] = np.nan
        vvf_temp = zfun.filt_hanning(vv, filtlen)
        Af[vname][:, ii] = vvf_temp
        vv_temp[ii] = np.nanmean(vvf_temp)
    Afm[vname] = vv_temp 
if doTS:
    TSf = dict()
    for vname in TS.keys():
        TSf[vname] = TS[vname].copy()
    TSfm = dict()    
    for vname in TS.keys():
        vv_temp = np.nan * np.ones(NP)
        for ii in range(NP):
            vv = TS[vname][:, ii]
            vv[iend_list[ii]:] = np.nan
            TSf[vname][:, ii] = zfun.filt_hanning(vv, filtlen)
            vv_temp[ii] = np.nanmean(vv)
        TSfm[vname] = vv_temp

# PREPARE FOR PLOTTING
# get the bathymetry
rpath = pmdir + 'roms/output/D2005_dia/'
fn = rpath + 'ocean_dia_1800.nc'
G, = zfun.get_basic_info(fn, getS=False, getT=False)
lon = G['lon_rho']
lat = G['lat_rho']
h = G['h']
# get the coastline
coast_file = pmdir + 'tools_data/geo_data/coast/pnw_coast_combined.mat'
import scipy.io
coast_mat = scipy.io.loadmat(coast_file)
coast_lon = coast_mat['lon'].squeeze()
coast_lat = coast_mat['lat'].squeeze()

# PLOTTING
plt.close()

fig = plt.figure(figsize=(16,8))

# map
ax = fig.add_subplot(121)
ax.contour(lon, lat, h, [2000, 1000, 500, 200, 100], colors='g')
ax.plot(coast_lon, coast_lat, '-k', linewidth=.5)
ax.axis([-126.5, -124, 47, 49.2])
zfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Winners & Losers')

# select a specific winner to explore
iiwin = iwin[0]

# plot track information

for ii in range(NP):
    if ii not in iwin:
        ax.plot(px[:, ii], py[:, ii], '-c', linewidth=.2)
        ax.plot(px[0, ii], py[0, ii], 'oc', markersize=3)
for ii in range(NP):
    if ii in iwin:
        ax.plot(px[:, ii], py[:, ii], '-b', linewidth=1)
        ax.plot(px[0, ii], py[0, ii], 'ob', markersize=10)
ax.plot(lon0, lat0,'oy', markersize=15)
ax.plot(px[:, iiwin], py[:, iiwin], '-r', linewidth=3)
ax.plot(px[0, iiwin], py[0, iiwin], '*r', markersize=20)

# time series

doCPG = True
if doCPG:    
    tag_list = ['t', 'cpg', 'visc']
    long_tag_list = ['Acceleration', 'Coriolis + PG', 'Friction']
    c_list = ['r', 'c', 'm']
    cpg_scale = 8
else:
    tag_list = ['t', 'cor', 'pg', 'visc']
    long_tag_list = ['Acceleration', 'Coriolis', 'Pressure Gradient', 'Friction']
    c_list = ['r', 'b', 'k', 'm']
    cpg_scale = 1

c_dict = dict(zip(tag_list, c_list))
lt_dict = dict(zip(tag_list, long_tag_list))

doFilt = True
if doFilt:
    ascale = 0.004 / cpg_scale
else:
    ascale = 0.02 / cpg_scale

ax1 = fig.add_subplot(222)
ax2 = fig.add_subplot(224)
for tag in tag_list:
    if doFilt:
        U = Af['U' + tag][:, iiwin]
        V = Af['V' + tag][:, iiwin]
    else:
        U = A['U' + tag][:, iiwin]
        V = A['V' + tag][:, iiwin]

    # rotate to point toward the mouth of JdF
    cal = np.cos(a0[iiwin])
    sal = np.sin(a0[iiwin])
    UU = U*cal + V*sal
    VV = V*cal - U*sal
    ax1.plot(td, UU, '-', color=c_dict[tag], linewidth=2, label=lt_dict[tag])
    ax2.plot(td, VV, '-', color=c_dict[tag], linewidth=2)
   
    if doFilt:    
        # add some mean vectors
        Ufm= Afm['U' + tag][iiwin]
        Vfm = Afm['V' + tag][iiwin]
        aa = ax.axis()
        print 'Warning - need to fix quiver call to handle scalar inputs: turn x into [x, x]'
        ax.quiver(aa[0] + .5, aa[2] + .5, Ufm, Vfm, color=c_dict[tag], units='y', scale=ascale, width=.02)
        
if doFilt:
    ax.text(aa[0] + .1, aa[2] + .2, 'PATH-AVERAGED FORCES')        
    
ax1.set_ylim(-ascale, ascale)
ax1.set_title('Acceleration along path toward goal (m s-2)')
ax1.grid()

if doFilt:
    ax1.text(.1, .1, 'Tidally Averaged', transform=ax1.transAxes)
else:
    ax1.text(.1, .1, 'Raw', transform=ax1.transAxes)
    
ax2.set_ylim(-ascale, ascale)
ax2.set_title('Acceleration across path')
ax2.grid()
ax2.set_xlabel('Time (days)')
    
ax1.legend(loc='upper left', fontsize=12)
        
plt.show()





