# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:03:12 2016

@author: PM5
"""

"""
Save results of Bradley Bartos' Rockfish experiments
to NetCDF.

"""

# setup
import os
import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import matfun
import pickle
import numpy as np
import pandas as pd
from datetime import datetime, timedelta

import netCDF4 as nc4

testing = True

if testing:
    sampsize = 100
    limit_days = True
else:
    sampsize = 1000 # 100000?
    limit_days = False

Ldir = Lfun.Lstart()

indir = '/data1/bbartos/LiveOcean_output/tracks/'
datadir = '/data1/bbartos/LiveOcean_data/tracker/'

# choose the run directory
# print('\n%s\n' % '** Choose mooring file to plot **')
# d_list_raw = os.listdir(indir)
# d_list = []
# for d in d_list_raw:
#     if 'MoSSea' in d:
#         d_list.append(d)
# d_list.sort()
# Ndt = len(d_list)
# for ndt in range(Ndt):
#     print(str(ndt) + ': ' + d_list[ndt])
# my_ndt = int(input('-- Input number -- '))
# dirname = d_list[my_ndt] + '/'

dirname = 'MoSSea_rockfish_rk4_ndiv1_forward_surfaceFalse_turbTrue_windage0_boundaryreflect/'

# create the list of run files
m_list_raw = os.listdir(indir + dirname)
m_list = []
for m in m_list_raw:
    m_list.append(m)
m_list.sort()
Npt = len(m_list)
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
my_ndt = int(input('-- Input number (99 for all) -- '))
if my_ndt == 99:
    pass
else:
    m_list = [m_list[my_ndt],]

# output directory
odir00 = Ldir['parent'] + 'ptools_output/'
Lfun.make_dir(odir00)
odir0 = odir00 + 'rockfish/'
Lfun.make_dir(odir0)
outdir = odir0
Lfun.make_dir(outdir)

# info for NetCDF output
name_unit_dict = {'lon':('Longitude','degrees'), 'lat':('Latitude','degrees'),
    'cs':('Fractional Z','Dimensionless'), 'ot':('Ocean Time','Seconds since 1/1/1970 UTC'),
    'z':('Z','m'), 'zeta':('Surface Z','m'), 'zbot':('Bottom Z','m'),
    'salt':('Salinity','Dimensionless'), 'temp':('Potential Temperature','Degrees C'),
    'u':('EW Velocity','meters s-1'), 'v':('NS Velocity','meters s-1'),
    'w':('Vertical Velocity','meters s-1'),
    'Uwind':('EW Wind Velocity','meters s-1'), 'Vwind':('NS Velocity','meters s-1'),
    'h':('Bottom Depth','m'), 'age':('Age','days')}

for inname in m_list:
    
    out_fn = outdir + inname + '.nc'
    try:
        os.remove(out_fn)
    except OSError:
        pass # assume error was because the file did not exist
    
    # compile list of day files
    p_list = os.listdir(indir + dirname + inname)
    p_list.sort()
    
    if limit_days == True:
        p_list = p_list[:5]
        
    # run through all days, concatenating the P dictionary in each
    counter = 0
    #P = dict()
    for p in p_list:
        if counter == 0:
            
            if not os.path.isfile(outdir + 'grid.nc'):
                # write a file of grid info
                dsh = nc4.Dataset('/pmr3/pmraid1/daves/runs/salish_2006_4/OUT/ocean_his_0025.nc')
                dsg = nc4.Dataset(outdir + 'grid.nc', 'w')
                # lists of variables to process
                dlist = ['xi_rho', 'eta_rho', 'xi_psi', 'eta_psi']
                vn_list2 = [ 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi', 'mask_rho', 'h']
                # Copy dimensions
                for dname, the_dim in dsh.dimensions.items():
                    if dname in dlist:
                        dsg.createDimension(dname, len(the_dim))
                # Copy variables
                for vn in vn_list2:
                    varin = dsh[vn]
                    vv = dsg.createVariable(vn, varin.dtype, varin.dimensions)
                    vv.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
                    vv[:] = dsh[vn][:]
                dsh.close()
                dsg.close()
                print('Wrote grid file.')
                sys.stdout.flush()
                
            # day 0 contains P, Ldir, and the grid data
            P, G, S, PLdir = pickle.load(open(indir + dirname + inname + '/' + p, 'rb'))
            NT, NP = P['lon'].shape
            ds = nc4.Dataset(out_fn, 'w')
            ds.createDimension('Time', None)
            ds.createDimension('Particle', NP)
            # Copy variables
            # for reference here is the full list
            vlist_full = ['lon', 'lat', 'cs', 'ot', 'z', 'zeta', 'zbot',
            'salt', 'temp', 'u', 'v', 'w', 'Uwind', 'Vwind', 'h', 'age']
            # but we will only carry these
            vlist = ['lon', 'lat', 'ot', 'z', 'h', 'age']
            for vn in vlist: #P.keys():
                if vn == 'ot':
                    vv = ds.createVariable(vn, float, ('Time'))
                else:
                    vv = ds.createVariable(vn, float, ('Time', 'Particle'))
                print(vn)
                print(P[vn].shape)
                if vn == 'age':
                    vv[:] = np.ones((NT,1)) * P[vn]
                else:
                    vv[:] = P[vn]
                print(ds[vn].shape)
                vv.long_name = name_unit_dict[vn][0]
                vv.units = name_unit_dict[vn]
            print('**')
            print(ds['lon'].shape)
            print('**')
            ds.close()
            
        else:
            # non-zero days only contain P and Ldir
            # first row overlaps with last row of previous day, so we remove it
            P, PLdir = pickle.load( open( indir + dirname + inname + '/' + p, 'rb' ) )
            # save the results
            ds = nc4.Dataset(out_fn, 'a')
            print(ds['lon'].shape)
            NTx, NPx = ds['lon'].shape
            print(NTx)
            for vn in vlist: #P.keys():
                if vn == 'ot':
                    ds[vn][NTx:] = P[vn][1:]
                else:
                    if vn == 'age':
                        ds[vn][NTx:,:] = np.ones((NT-1,1)) * P[vn]
                    else:
                        ds[vn][NTx:,:] = P[vn][1:,:]
                print(vn)
                print(ds[vn].shape)
            ds.close()
        counter += 1
        #print('Finished ' + p)
        sys.stdout.flush()
