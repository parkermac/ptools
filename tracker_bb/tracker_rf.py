"""
Code for particle tracking of rockfish larvae.
Derivative of tracker_b.py

With rockfish experiment, larvae are released at one location 
over two months and following the lunar cycle. 
Then, experiment continues to run for four months.
Will run for years listed in yr_list, default is 2006.
If changing the years, update the lun_dt_ind as well.

Designed to read csv files. One needs to contain the columns:
"Experiment", "model", "species", "lat", "lon", and "depth (m)",.
The other must contain: "days", and "# of particles".

The experiment will be chosen by either an argument >tracker_rf.py -ex 1_1
or a user input after script has started to run

Post-Processing Alterations:
To add daily mortality, use mort_daily.py
To impose time limit (6 months) on particles after release, use mort_age.py
"""

#%% setup
import argparse
import numpy as np
import pandas as pd
import netCDF4 as nc
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
import trackfun
reload(trackfun)

import warnings
warnings.simplefilter('error')

# some run specifications
ic_name = 'rockfish'
dir_tag = 'forward' # 'forward' or 'reverse'
method = 'rk4' # 'rk2' or 'rk4'
surface = False # Boolean, True for trap to surface
turb = True # Boolean, True to include vertical turbulence
windage = 0 # a small number >= 0
ndiv = 1 # number of divisions to make between saves for the integration
        # e.g. if ndiv = 3 and we have hourly saves, we use a 20 minute step
        # for the integration (but still only report fields hourly)
bound = 'reflect' # either 'stop' for halting at coast or 'reflect to bounce off

# create possible inputs
parser = argparse.ArgumentParser()
parser.add_argument('-e', '--experiment', nargs='?', type=str)
parser.add_argument('-t', '--testing', action='store_true') # Use to limit the number of years, particles, and timesteps, default is False
args = parser.parse_args()
testing = args.testing

# retrieve experimental data
exdf = pd.read_csv(Ldir['data'] + 'tracker/rockfish_latlon.csv', index_col = 0)
pardf = pd.read_csv(Ldir['data'] + 'tracker/rockfish_partruition.csv', index_col=0)

# choose experiment if not already provided as an argument
if args.experiment != None:
    exrow = args.experiment
else:
    print('\n%s\n' % '** Choose Experiment **')
    for ind in exdf.index:
        print(ind+'  '+str(exdf.loc[ind][0])+'  '+str(exdf.loc[ind][1]))
    exrow = str(input('-- Input Experiment -- '))

# set number of particles
if testing:
    NP0 = 100
else:
    NP0 = exdf.loc[exrow]['# of particles']

# set particle initial locations
plon00 = np.array([exdf.loc[exrow]['lon']])
plat00 = np.array([exdf.loc[exrow]['lat']])
dep_orig = exdf.loc[exrow]['depth (m)']
# depth array has full number of larvae
dep00 = np.linspace(dep_orig, dep_orig, NP0)

# set model and species
gtagex = exdf.loc[exrow]['model']
species = exdf.loc[exrow]['species']

# set number of days
if testing:
    days_to_track = 5
else:
    days_to_track = 180 # standard run will be 6 months

# years to track
yr_list = [2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009]
# start dates for each year coincide with new moon
if species == 'canary':
    lunar_dt_dict = dict(zip(yr_list, pd.to_datetime(['2002/01/13', 
                  '2003/01/02', '2004/01/21', '2005/01/10', '2006/01/29',
                  '2007/01/18', '2008/01/08', '2009/01/25'])))
elif species == 'yelloweye':
    lunar_dt_dict = dict(zip(yr_list, pd.to_datetime(['2002/05/12', 
                  '2003/05/01', '2004/05/18', '2005/05/08', '2006/05/26', 
                  '2007/05/16', '2008/05/05', '2009/05/24'])))

# choose year for testing or full list
if testing:
    yr_list = [2006,]
else:
    pass

# MoSSea model only has one year
if gtagex == 'MoSSea':
    yr_list = [2006,]

# make the full IC vectors, which will have equal length
NSP = len(dep00)
NXYP = len(plon00)
plon0 = plon00.reshape(NXYP,1) * np.ones((NXYP,NSP))
plat0 = plat00.reshape(NXYP,1) * np.ones((NXYP,NSP))
dep0 = dep00.reshape(1,NSP) * np.ones((NXYP,NSP))
plon0 = plon0.flatten()
plat0 = plat0.flatten()
dep0 = dep0.flatten()

# create initial depth ranges
dep_min = np.ones(len(plon0)) * -20
dep_max = np.ones(len(plon0)) * -100
dep_range = (dep_min, dep_max)

# create array for tracking particle age
age = np.zeros(len(plon0))

# save some things in Ldir
Ldir['gtagex'] = gtagex
Ldir['ic_name'] = ic_name
Ldir['dir_tag'] = dir_tag
Ldir['method'] = method
Ldir['surface'] = surface
Ldir['turb'] = turb
Ldir['windage'] = windage
Ldir['ndiv'] = ndiv
Ldir['experiment'] = exrow
Ldir['days_to_track'] = days_to_track

# make sure the output directory exists
outdir0 = Ldir['LOo'] + 'tracks/'
Lfun.make_dir(outdir0)

#%% Tracking

# track in year loop
for yr in yr_list:
    print('Working on ' + str(yr))
    print(' - Starting on ' + time.asctime())
    sys.stdout.flush()

    idt0 = lunar_dt_dict[yr]
    
    # make the output directory to store multiple experiments
    outdir = (outdir0 +
        Ldir['gtagex'] + '_' + Ldir['ic_name'] + '_' + Ldir['method'] +
        '_' + 'ndiv' + str(Ldir['ndiv']) + '_' + Ldir['dir_tag'] +
        '_' + 'surface' + str(Ldir['surface']) + '_' + 'turb' + 
        str(Ldir['turb']) + '_' + 'windage' + str(Ldir['windage']) + 
        '_boundary' + bound + '_nodepthchange/')
    Lfun.make_dir(outdir)

    # split the calculation up into one-day chunks    
    for nd in range(Ldir['days_to_track']):
        
        idt = idt0 + timedelta(days=nd)
        
        # make sure out file list starts at the start of the day
        if nd > 0: 
            fn_last = fn_list[-1]
        fn_list = trackfun.get_fn_list(idt, Ldir, yr=yr)
        if nd > 0:
            fn_list = [fn_last] + fn_list
        
        # set dates
        T0 = zrfun.get_basic_info(fn_list[0], only_T=True)
        Tend = zrfun.get_basic_info(fn_list[-1], only_T=True)
        Ldir['date_string0'] = datetime.strftime(T0['tm'],'%Y.%m.%d')
        Ldir['date_string1'] = datetime.strftime(Tend['tm'],'%Y.%m.%d')
        
        # change depth ranges based on age
        juv = age >= 40
        juv_min = np.ones(len(plon0)) * -50
        juv_max = np.ones(len(plon0)) * -100
#        if sum(juv) != 0:
#            dep_range[0][juv] = juv_min[juv] 
#            dep_range[1][juv] = juv_max[juv]
        
        if nd == 0: # first day
        
            # make the output subdirectory for this year
            outdir1 = (outdir + Ldir['ic_name'] + '_' + str(yr) + '_Experiment_' + exrow + '/')
            Lfun.make_dir(outdir1, clean=True)
            
            # checking initial locations
            # pull start dataset
            ds0 = nc.Dataset(fn_list[0])
            
            # make vectors to feed to interpolant maker
            G = zrfun.get_basic_info(fn_list[0], only_G=True)
            S = zrfun.get_basic_info(fn_list[0], only_S=True)
            R = dict()
            R['rlonr'] = G['lon_rho'][0,:].squeeze()
            R['rlatr'] = G['lat_rho'][:,0].squeeze()
            R['rlonu'] = G['lon_u'][0,:].squeeze()
            R['rlatu'] = G['lat_u'][:,0].squeeze()
            R['rlonv'] = G['lon_v'][0,:].squeeze()
            R['rlatv'] = G['lat_v'][:,0].squeeze()
            R['rcsr'] = S['Cs_r'][:]
            R['rcsw'] = S['Cs_w'][:]
            
            # when initial position is on the boundary, need to move one 
            # grid cell towards deeper area so velocities are not zero
            pcs_temp = np.array([0])
            plon0, plat0 = trackfun.change_position(ds0, plon0, plat0, pcs_temp, R, surface)
        
            # create initial pcs from depths
            ZH0 = trackfun.get_V(['zeta', 'h'], ds0, plon0, plat0, pcs_temp, R, surface)
            Tot_Dep = ZH0[0,0] + ZH0[0,1]
            pcs0 = -dep0/Tot_Dep
            
            # do the tracking
            tt0 = time.time()
            P, G, S = trackfun.get_tracks(fn_list, plon0, plat0, pcs0, dir_tag,
                                  method, surface, turb, ndiv, windage, 
                                  bound=bound, dep_range=dep_range)
            
            # reset locations of particles not released
            if testing:
                pass
            else:
                npart = pardf['running sum'].iloc[nd]
                P['lon'][:,npart:] = plon0[npart:]
                P['lat'][:,npart:] = plat0[npart:]
                P['cs'][:,npart:] = pcs0[npart:]
            P['age'] = age
                                  
            # save the results
            outname = 'day_' + ('00000' + str(nd))[-5:] + '.p'
            pickle.dump( (P, G, S, Ldir) , open( outdir1 + outname, 'wb' ) )
            print(' - Took %0.1f sec to produce ' %(time.time()-tt0) + outname)
            sys.stdout.flush()   
            
        else: # subsequent days
            tt0 = time.time()
            # get initial condition
            plon0 = P['lon'][-1,:]
            plat0 = P['lat'][-1,:]
            pcs0 = P['cs'][-1,:]
            # do the tracking
            P, G, S = trackfun.get_tracks(fn_list, plon0, plat0, pcs0, dir_tag,
                                  method, surface, turb, ndiv, windage, 
                                  bound=bound, dep_range=dep_range)
            
            # reset locations of particles not released
            if testing:
                pass
            else:
                if nd < 59:
                    npart = pardf['running sum'].iloc[nd]
                    P['lon'][:,npart:] = plon0[npart:]
                    P['lat'][:,npart:] = plat0[npart:]
                    P['cs'][:,npart:] = pcs0[npart:]
                    # update age array with running sum of particles
                    age[:npart] = age[:npart] + 1
                else:
                    age = age + 1
            P['age'] = age
            
            # save the results
            outname = 'day_' + ('00000' + str(nd))[-5:] + '.p'
            pickle.dump( (P, Ldir) , open( outdir1 + outname, 'wb' ) )
            print(' - Took %0.1f sec to produce ' %(time.time()-tt0) + outname)
            sys.stdout.flush() 

    print('Results saved to:\n' + outdir)
    print(50*'*')
    sys.stdout.flush() 
