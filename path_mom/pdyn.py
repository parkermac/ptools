"""
Code to analyze the dynamics of flow along particle trajectories.

This is the first piece of code to run (after ROMS and particulator). It reads
in a particulator file - assumed to consist of a single release of a bunch of
particles.  Then for each location of each particle it goes into the ROMS
diagnostic output files and interpolates in time and space to get all the
diagnostics.  These are added with their native names to the dict "p" - into
which we have already read the particulator data, and then saved using pickle.
"""

# USER: set values
pmdir = '/Users/PM3/Documents/'
pfn = pmdir + 'tools_output/particulator_out/jdf_inflow/yd.75.nc' 
coast_file = pmdir + 'tools_data/geo_data/coast/pnw_coast_combined.mat'

testing = False

# IMPORTS
#import matplotlib.pyplot as plt
import numpy as np
import zfun; reload(zfun)

# make a list of roms files
rpath = '../../roms/output/D2005_dia/'
nn_list = range(1800, 1923 + 1) # user responsible to see that these exist!
rfn_list = []
for nn in nn_list:
    ns = ('0000' + str(nn))[-4:]
    fn = rpath + 'ocean_dia_' + ns + '.nc'
    rfn_list.append(fn)

# get the roms grid info
G, S = zfun.get_basic_info(rfn_list[0], getT=False)

# get the time    
rNT = len(nn_list)
rt = np.zeros((rNT,1))
ii = 0
for rfn in rfn_list:
    [T] = zfun.get_basic_info(rfn, getG=False, getS=False)
    rt[ii] = float(T['ocean_time'])
    ii += 1
rtday = rt/86400.

# get the particulator data
import netCDF4 as nc
ds = nc.Dataset(pfn)
p = dict()
for vn in ds.variables:
    # print ds.variables[vn] # to see info
    p[vn] =  ds.variables[vn][:]
ds.close()

pNT, NP = p['t'][:].shape
pt = p['t'][:,0] # time in seconds from start of year (like tsec)

if testing:
    pNT = 5

itp_list = zfun.get_interpolant(pt, rt)
# now each element in itp_list is a three-element "interpolant" tuple, with:
# index into rfn_list of history file at time before
# index into rfn_list of history file at time after
# fraction of the way between the two

rlonr = G['lon_rho'][0,:].squeeze()
rlatr = G['lat_rho'][:,0].squeeze()
rlonu = G['lon_u'][0,:].squeeze()
rlatu = G['lat_u'][:,0].squeeze()
rlonv = G['lon_v'][0,:].squeeze()
rlatv = G['lat_v'][:,0].squeeze()
rcs = S['Cs_r'][:]

ts_list = ['rate','xadv','yadv','vadv','hdiff','vdiff']
uv_list = ['accel','xadv','yadv','vadv','cor','prsgrd','vvisc']
for ts in ts_list[:]:
    p['salt_' + ts] = np.zeros((pNT, NP))
    p['temp_' + ts] = np.zeros((pNT, NP))
for uv in uv_list[:]:
    p['u_' + uv] = np.zeros((pNT, NP))
    p['v_' + uv] = np.zeros((pNT, NP))

import time 
# start a timer
tt0 = time.time()  

tt = 0
for itp in itp_list[:pNT]: # use [:5] e.g. to speed up testing
    print ' ** working on time ' + str(tt) + ' out of ' + str(pNT)
    
    fn0 = rfn_list[itp[0]]
    fn1 = rfn_list[itp[1]]
   
    plon = p['lon'][tt, :]
    plat = p['lat'][tt, :]
    pcs = p['cs'][tt, :]
    
    ilonr_list = zfun.get_interpolant(plon, rlonr)
    ilatr_list = zfun.get_interpolant(plat, rlatr)
    
    ilonu_list = zfun.get_interpolant(plon, rlonu)
    ilatu_list = zfun.get_interpolant(plat, rlatu)
    
    ilonv_list = zfun.get_interpolant(plon, rlonv)
    ilatv_list = zfun.get_interpolant(plat, rlatv)
    
    ics_list = zfun.get_interpolant(pcs, rcs)
           
    ds0 = nc.Dataset(fn0)
    ds1 = nc.Dataset(fn1)
    
    for pp in range(NP):
        
        ilonr = ilonr_list[pp]
        ilatr = ilatr_list[pp]
        
        ilonu = ilonu_list[pp]
        ilatu = ilatu_list[pp]
        
        ilonv = ilonv_list[pp]
        ilatv = ilatv_list[pp]
        
        ics = ics_list[pp]
        if ics[1] == S['N']: # sometimes ics gets out of bounds by 1
            iics = tuple([ics[0]-1, ics[1]-1, ics[2]])
        else:
            iics = ics
            
        for ts in ts_list[:]:
                        
            varname = 'salt_' + ts            
            ival = zfun.interpolate4D(ds0, ds1, varname, itp, iics, ilatr, ilonr)           
            p[varname][tt, pp] = ival
            
            varname = 'temp_' + ts            
            ival = zfun.interpolate4D(ds0, ds1, varname, itp, iics, ilatr, ilonr)           
            p[varname][tt, pp] = ival           
            
        for uv in uv_list[:]:
            
            varname = 'u_' + uv            
            ival = zfun.interpolate4D(ds0, ds1, varname, itp, iics, ilatu, ilonu)           
            p[varname][tt, pp] = ival
            
            varname = 'v_' + uv            
            ival = zfun.interpolate4D(ds0, ds1, varname, itp, iics, ilatv, ilonv)           
            p[varname][tt, pp] = ival
    
    ds0.close()
    ds1.close()
    
    tt += 1
    
print str(time.time() - tt0) + ' sec to get ' + str(tt) + ' rates'
# RESULT takes 385 sec to get all temp and salt rates (361 times, 105 particles)
# or 6.5 minutes per 5-day experiment - a full year would take overnight

import cPickle
pout = open(pmdir + 'tools_output/pydev_out/p75.pkl', 'wb')
cPickle.dump(p, pout)
pout.close()
# RESULT this file is 30MB for this experiment






