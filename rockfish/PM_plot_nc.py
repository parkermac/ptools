# -*- coding: utf-8 -*-
"""
Plot results of tracker using the new .nc files.
"""

# setup
import os
import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

plp = os.path.abspath('../../LiveOcean/plotting')
if plp not in sys.path:
    sys.path.append(plp)
import pfun

import numpy as np
import netCDF4 as nc
import pandas as pd

from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

if Ldir['lo_env'] == 'pm_mac': # mac version
    pass
elif Ldir['lo_env'] == 'pm_fjord': # fjord version
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt

# create the list of run files
indir = odir00 = Ldir['parent'] + 'ptools_output/rockfish/'
ex_list_raw = os.listdir(indir)
ex_list = []
for ex in ex_list_raw:
    if ('.nc' in ex) and ('grid' not in ex):
        ex_list.append(ex)
ex_list.sort()

plotdir = indir + 'plots/'
Lfun.make_dir(plotdir)
save_plots = True

# testing
ex_list = ex_list[-12:]
#ex = 'rockfish_2006_Experiment_4_4.nc'

# retrieve experimental data
if Ldir['lo_env'] == 'pm_fjord':
    datadir = '/data1/bbartos/LiveOcean_data/tracker/'
elif Ldir['lo_env'] == 'pm_mac':
    datadir = Ldir['parent'] + 'ptools_data/rockfish/'
exdf = pd.read_csv(datadir + 'rockfish_latlon.csv', index_col = 0)

plt.close('all')

for ex in ex_list:
    
    try:

        print('Plotting ' + ex)
        # get data
        ds = nc.Dataset(indir + ex)
        # tracks are stored (time, particle)
        Lon = ds['lon'][:]
        Lat = ds['lat'][:]
        Z = ds['z'][:]
        H = ds['h'][:]
        Age = ds['age'][:]
        ot = ds['ot'][:]
        ds.close()
        
        # subsample the number of particles
        pstep = 100 # use 100 typically
        lon = Lon[:,::pstep]
        lat = Lat[:,::pstep]
        z = Z[:,::pstep]
        h = H[:,::pstep]
        age = Age[:,::pstep]
        NT, NP = lon.shape
        
        # find time indices that define an
        # age range (e.g. 0-120 days)
        ndays0 = 0
        ndays1 = 120
        age_ind0 = (age == ndays0).sum(axis=0)
        age_ind1 = (age <= ndays1).sum(axis=0)
        
        # repackage to all have the same age
        nt = age_ind1[0] - age_ind0[0]
        llon = np.nan * np.ones((nt, NP))
        llat = np.nan * np.ones((nt, NP))
        hh  = np.nan * np.ones((nt, NP))
        zz  = np.nan * np.ones((nt, NP))
        for npt in range(NP):
            nt0 = age_ind0[npt]
            nt1 = nt0 + nt
            llon[:,npt] = lon[nt0:nt1, npt]
            llat[:,npt] = lat[nt0:nt1, npt]
            hh[:,npt] = h[nt0:nt1, npt]
            zz[:,npt] = z[nt0:nt1, npt]
        
        # Drop particles that end up deeper than 50 m
        # because the shallow ones are typically stuck on land.
        # Also drop particles that escape the domain.
        mask1 = hh[-1,:] > 50
        mask2 = llon.min(axis=0) > -126.5
        mask3 = llat.min(axis=0) > 45.5
        mask = mask1 & mask2 & mask3
        llon = llon[:,mask]
        llat = llat[:,mask]
        hh = hh[:,mask]
        zz = zz[:,mask]
        
        # get strings that describe the experiment, species, and release location
        exx = ex.split('_')
        ex1 = exx[-2]
        ex2 = exx[-1]
        ex2 = ex2.replace('.nc','')
        exind = ex1 + '_' + ex2
        species = exdf.loc[exind, 'species']
        location = exdf.loc[exind, 'Site']
        
        # PLOTTING
        fig = plt.figure(figsize=(13,7))
        
        # MAP
        ax = fig.add_subplot(1,2,1)
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        pfun.add_coast(ax)
        # automatically plot region of particles, with padding
        pad = .02
        aa = [llon.min() - pad, llon.max() + pad,
            llat.min() - pad, llat.max() + pad]
        #aa = [-123.5, -122, 47, 48.5]
        ax.axis(aa)
        pfun.dar(ax)
        ax.grid()
        ax.plot(llon, llat, '-', linewidth=0.5, alpha=0.5)
        ax.plot(llon[0,0], llat[0,0], 'or', markersize=10) # start
        ax.plot(llon[-1,:],llat[-1,:], 'ob', markersize=3) # end
        ax.set_title('%s : %s : %s' % (exind, species.title(), location))
        
        # TIME SERIES
        age_days = np.linspace(0, nt/24, nt)
        ax = fig.add_subplot(2,2,2)
        ax.plot(age_days, zz,'-', alpha=0.25)
        ax.set_ylabel('Z (m)')
        ax.set_xlabel('Age (days)')
        
        # save or plot figures
        if save_plots:
            out_fn = plotdir + ex.strip('.nc') + '.png'
            plt.savefig(out_fn)
            plt.close()
        else:
            plt.show()
        
    except:
        # Experiment 4_15 doesn't plot right (no particles)
        print('-- could not plot ' + ex)
        pass


   

