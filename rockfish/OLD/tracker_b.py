"""
Code for particle tracking.  Meant to be fast and simple, but capable
of tracking backward in time.

Designed for ROMS output with plaid lat, lon grids.

Recoded the integration scheme around 7/15/2016 to use new version of
get_interpolant.py, and to better handle many particles.

PERFORMANCE: With the new fast version is took about 12 seconds
for a day of integration for 6-600 particles, and only increased to
19 seconds for 6000 particles.  The old slow version was faster
for < 100 particles, but otherwise became very slow, scaling linearly
with the number of particles.  These tests were all with the cascadia1 grid.

The design philosophy is that this should be capable of handling both
LiveOcean and older versions, like from PNWTOX, of how files are stored.
It should also work on different platforms (my mac and fjord at this time).

This program is mainly a DRIVER where you supply:
    - which run to work on
    - particle initial locations
    - a list of starting times
    - duration in days to track
    - forward or backward in time
    - some other flags

"""

#%% setup
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import time
import pickle

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()
# Note that we override parts of the Ldir logic,
# using Ldir['gtagex'] to identify a run.

from importlib import reload
import zrfun 
reload(zrfun)
# using trackfun_b here, if fully integrated with trackfun:
# just delete this comment and simplify the statement
import trackfun_b as trackfun
reload(trackfun)

#%% ************ USER INPUT **************************************

# some run specifications
gtagex = 'cascadia1_base_lobio1' # 'cascadia1_base_lobio1', 'D2005_his','C2009'
ic_name = 'turbdemo' # 'jdf', 'cr', 'deadbirds', 'akashiwo', 'test', 'turbdemo'
dir_tag = 'forward' # 'forward' or 'reverse'
method = 'rk4' # 'rk2' or 'rk4'
surface = False # Boolean, True for trap to surface
turb = True # Boolean, True to include vertical turbulence
stagger = True # staggering particles with fn_mask
windage = 0 # a small number >= 0
ndiv = 1 # number of divisions to make between saves for the integration
        # e.g. if ndiv = 3 and we have hourly saves, we use a 20 minute step
        # for the integration (but still only report fields hourly)

# set run time information (multiple experiments)
# always start on a day (no hours)
if Ldir['parent'] == '/Users/PM5/Documents/':
    # mac version
    if gtagex == 'cascadia1_base_lo1':
        dt_first_day = datetime(2015,9,19)
        number_of_start_days = 1
        days_between_starts = 1
        days_to_track = 1
    elif gtagex == 'cascadia1_base_lobio1':
        dt_first_day = datetime(2015,9,16)
        number_of_start_days = 1
        days_between_starts = 1
        days_to_track = 1
    elif gtagex == 'D2005_his':
        dt_first_day = datetime(2005,3,17)
        number_of_start_days = 3
        days_between_starts = 1
        days_to_track = 2
    elif gtagex == 'C2009':
        dt_first_day = datetime(2009,8,1)
        number_of_start_days = 4
        days_between_starts = 7
        days_to_track = 7

elif Ldir['parent'] == '/data1/parker/':
    # fjord version
    if gtagex == 'cascadia1_base_lo1':
        dt_first_day = datetime(2014,11,1)
        number_of_start_days = 48
        days_between_starts = 3
        days_to_track = 7        
    elif gtagex == 'C2009':
        dt_first_day = datetime(2009,8,1)
        number_of_start_days = 4
        days_between_starts = 7
        days_to_track = 7

elif Ldir['parent'] == 'C:/Users/Bradley/Documents/Research Work/Parker/':
    # Bradley's PC
    if gtagex == 'cascadia1_base_lobio1': 
        dt_first_day = datetime(2016,1,1)
        number_of_start_days = 1
        days_between_starts = 2
        days_to_track = 5
    elif gtagex == 'C2009':
        dt_first_day = datetime(2009,8,24)
        number_of_start_days = 4
        days_between_starts = 7
        days_to_track = 7

elif Ldir['parent'] == '/home/bbartos/':
    # Bradley's fjord
    if gtagex == 'cascadia1_base_lobio1':
        dt_first_day = datetime(2016,7,23)
        number_of_start_days = 1
        days_between_starts = 3
        days_to_track = 5
    elif gtagex == 'C2009':
        dt_first_day = datetime(2009,8,24)
        number_of_start_days = 4
        days_between_starts = 7
        days_to_track = 7
                       
# set particle initial locations, all numpy arrays
#
# first create three vectors of initial locations
# plat00 and plon00 should be the same length,
# and the length of pcs00 is however many vertical positions you have at
# each lat, lon

if ic_name == 'jdf':
    plon00 = np.array([-124.65])
    plat00 = np.array([48.48])
    pcs00 = np.linspace(-.95,-.05,20)
elif ic_name == 'cr':
    plon00 = np.array([-123.9])
    plat00 = np.array([46.22])
    pcs00 = np.linspace(-.95,-.05,20)
elif ic_name in ['deadBirds', 'test']:
    lonvec = np.linspace(-127, -123.9, 20)
    latvec = np.linspace(43.5, 49.5, 30)
    lonmat, latmat = np.meshgrid(lonvec, latvec)
    plon00 = lonmat.flatten()
    plat00 = latmat.flatten()
    pcs00 = np.array([-.05])
elif ic_name == 'akashiwo':
    lldf = pd.read_csv(Ldir['data']+'tracker/'+ic_name+'.csv', index_col=0)
    plon00 = lldf['lon']
    plat00 = lldf['lat']
    pcs00 = np.array([-0.05])
elif ic_name == 'turbdemo':
    plon00 = np.array([-124.4])
    plat00 = np.array([47.25])
    pcs00 = np.linspace(-.05, -.05, 10)

if len(plon00) != len(plat00):
    print('Problem with length of initial lat, lon vectors')
    sys.exit()
# ********* END USER INPUT *************************************


# save some things in Ldir
Ldir['gtagex'] = gtagex
Ldir['ic_name'] = ic_name
Ldir['dir_tag'] = dir_tag
Ldir['method'] = method
Ldir['surface'] = surface
Ldir['turb'] = turb
Ldir['windage'] = windage
Ldir['ndiv'] = ndiv
Ldir['days_to_track'] = days_to_track

# make the full IC vectors, which will have equal length
# (one value for each particle)
NSP = len(pcs00)
NXYP = len(plon00)
plon0 = plon00.reshape(NXYP,1) * np.ones((NXYP,NSP))
plat0 = plat00.reshape(NXYP,1) * np.ones((NXYP,NSP))
pcs0 = pcs00.reshape(1,NSP) * np.ones((NXYP,NSP))
plon0 = plon0.flatten()
plat0 = plat0.flatten()
pcs0 = pcs0.flatten()

# make sure the output directory exists
outdir0 = Ldir['LOo'] + 'tracks/'
Lfun.make_dir(outdir0)

# make the list of start days (datetimes)
idt_list = []
dt = dt_first_day
for nic in range(number_of_start_days):
    idt_list.append(dt)
    dt = dt + timedelta(days_between_starts)

#%% step through the experiments (one for each start day)
for idt0 in idt_list:

    # make the output directory to store multiple experiments
    outdir = (outdir0 +
        Ldir['gtagex'] + '_' + Ldir['ic_name'] + '_' + Ldir['method'] +
        '_' + 'ndiv' + str(Ldir['ndiv']) + '_' + Ldir['dir_tag'] +
        '_' + 'surface' + str(Ldir['surface']) + '_' + 'turb' + 
        str(Ldir['turb']) + '_' + 'windage' + str(Ldir['windage']) + 
        '_' + datetime.strftime(idt0,'%Y.%m.%d') + '_' + 
        str(Ldir['days_to_track']) + 'days/')
    Lfun.make_dir(outdir)

    # split the calculation up into one-day chunks    
    for nd in range(Ldir['days_to_track']):
        
        idt = idt0 + timedelta(days=nd)
        
        # make sure out file list starts at the start of the day
        if nd > 0: 
            fn_last = fn_list[-1]
        fn_list = trackfun.get_fn_list(idt, Ldir)
        if nd > 0:
            fn_list = [fn_last] + fn_list
        
        # get starting and ending dates
        T0 = zrfun.get_basic_info(fn_list[0], only_T=True)
        Tend = zrfun.get_basic_info(fn_list[-1], only_T=True)
        Ldir['date_string0'] = datetime.strftime(T0['tm'],'%Y.%m.%d')
        Ldir['date_string1'] = datetime.strftime(Tend['tm'],'%Y.%m.%d')
        
        if nd == 0: # first day
            print(50*'*')
            print('Calculating tracks from ' + Ldir['date_string0'])
            
            # remove particles on land for gridded experiments
            if ic_name in ['deadBirds', 'test']:
                # get basic info
                G = zrfun.get_basic_info(fn_list[0], only_G=True)
                S = zrfun.get_basic_info(fn_list[0], only_S=True)
                ds0 = nc4.Dataset(fn_list[0])
            
                # make vectors to feed to interpolant maker
                R = dict()
                R['rlonr'] = G['lon_rho'][0,:].squeeze()
                R['rlatr'] = G['lat_rho'][:,0].squeeze()
                R['rlonu'] = G['lon_u'][0,:].squeeze()
                R['rlatu'] = G['lat_u'][:,0].squeeze()
                R['rlonv'] = G['lon_v'][0,:].squeeze()
                R['rlatv'] = G['lat_v'][:,0].squeeze()
                R['rcsr'] = S['Cs_r'][:]
                R['rcsw'] = S['Cs_w'][:]
                
                # pull salinity field with nan on land
                SALT = trackfun.get_V(['salt'], ds0, plon0, plat0, pcs0, R, surface)
                SALT = SALT.flatten()
                plon0 = plon0[~np.isnan(SALT)]
                plat0 = plat0[~np.isnan(SALT)]
                pcs0 = pcs0[~np.isnan(SALT)]
            
            # make the output subdirectory for this experiment
            outdir1 = (outdir + Ldir['ic_name'] + '_' + Ldir['date_string0'] + 
                        '_' + str(Ldir['days_to_track']) + 'days/')
            Lfun.make_dir(outdir1, clean=True)
            
            # do the tracking
            tt0 = time.time()
            P, G, S = trackfun.get_tracks(fn_list, plon0, plat0, pcs0, dir_tag,
                                          method, surface, turb, ndiv, windage)
            
            # save the results
            outname = 'day_' + ('00000' + str(nd))[-5:] + '.p'
            pickle.dump( (P, G, S, Ldir) , open( outdir1 + outname, 'wb' ) )
            print(' - Took %0.1f sec to produce ' %(time.time()-tt0) + outname)
            
        else: # subsequent days
            tt0 = time.time()
            # get initial condition from previous day
            plon0 = P['lon'][-1,:]
            plat0 = P['lat'][-1,:]
            pcs0 = P['cs'][-1,:]
            
            # do the tracking
            P, G, S = trackfun.get_tracks(fn_list, plon0, plat0, pcs0, dir_tag,
                                          method, surface, turb, ndiv, windage)
            
            # save the results
            outname = 'day_' + ('00000' + str(nd))[-5:] + '.p'
            pickle.dump( (P, Ldir) , open( outdir1 + outname, 'wb' ) )
            print(' - Took %0.1f sec to produce ' %(time.time()-tt0) + outname) 

    print('Results saved to:\n' + outdir)
    print(50*'*')
