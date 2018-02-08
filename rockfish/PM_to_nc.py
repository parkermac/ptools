# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:03:12 2016

@author: PM5
"""

"""
Save results of Bradley Bartos' Rockfish experiments
to NetCDF.  Can limit the variables to save space.

Takes about 2 minutes to make the .nc file for an experiment (all 180 days)
and the resulting file is 17GB, which takes me 1/2 hour to download at home.

"""

# setup
import os
import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import pickle
import numpy as np
import time
import argparse
import netCDF4 as nc4


# optional command line arguments, can be input in any order
parser = argparse.ArgumentParser()
parser.add_argument('-n', '--exp_num', nargs='?', type=int, default=0)
args = parser.parse_args()

Ldir = Lfun.Lstart()

daylim = 6
if Ldir['lo_env'] == 'pm_fjord':
    limit_days = False
    pstep = 10 # subsampling factor for number of particles
    # here we hardwire a single SET OF EXPERIMENTS because I believe these are the only
    # ones we want to work with (all other things like mortality can be
    # treated in post-processing of these NetCDF files)
    expdir0 = ('/data1/bbartos/LiveOcean_output/tracks/' +
        'MoSSea_rockfish_rk4_ndiv1_forward_surfaceFalse_turbTrue_windage0_boundaryreflect/')
    # create the list of experiment directories
    exp_list_raw = os.listdir(expdir0)
    exp_list = []
    for m in exp_list_raw:
        if 'Experiment' in m:
            exp_list.append(m)
    exp_list.sort()
elif Ldir['lo_env'] == 'pm_mac':
    limit_days = True
    pstep = 100
    expdir0 = Ldir['parent'] + 'ptools_data/rockfish/'
    exp_list = ['rockfish_2006_Experiment_3_1']
    args.exp_num = 0

if args.exp_num == 99:
    # leaves expdir_list intact so it does the all
    pass
else:
    # just a single item in the list
    exp_list = [exp_list[args.exp_num]]
    
for m in exp_list:
    print(m)
print('')

# output directory
odir0 = Ldir['parent'] + 'ptools_output/'
Lfun.make_dir(odir0)
outdir = odir0 + 'rockfish/'
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

for exp in exp_list:
    print('** Working on ' + exp)
    tt0 = time.time()
    out_fn = outdir + exp + '.nc'
    print(' ' + out_fn)
    try:
        os.remove(out_fn)
    except OSError:
        pass # assume error was because the file did not exist
    # compile list of day files
    p_list = os.listdir(expdir0 + exp)
    p_list.sort()
    if limit_days == True:
        p_list = p_list[:daylim]
    counter = 0
    for p in p_list:
        print(p)
        sys.stdout.flush()
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
            P, G, S, PLdir = pickle.load(open(expdir0 + exp + '/' + p, 'rb'))
            NT, NP = P['lon'][:,::pstep].shape
            ds = nc4.Dataset(out_fn, 'w')
            ds.createDimension('Time', None)
            ds.createDimension('Particle', NP)
            # Copy variables
            # for reference here is the full list
            vlist_full = ['lon', 'lat', 'cs', 'ot', 'z', 'zeta', 'zbot',
            'salt', 'temp', 'u', 'v', 'w', 'Uwind', 'Vwind', 'h', 'age']
            # but we will only carry these
            vlist = ['lon', 'lat', 'ot', 'z', 'h', 'age']
            for vn in vlist:
                if vn == 'ot':
                    vv = ds.createVariable(vn, float, ('Time'))
                else:
                    vv = ds.createVariable(vn, float, ('Time', 'Particle'))
                if vn == 'ot':
                    vv[:] = P[vn]
                elif vn == 'age':
                    vv[:] = np.ones((NT,1)) * P[vn][::pstep]
                else:
                    vv[:] = P[vn][:,::pstep]
                vv.long_name = name_unit_dict[vn][0]
                vv.units = name_unit_dict[vn][1]
            ds.close()
        else:
            # non-zero days only contain P and Ldir
            P, PLdir = pickle.load( open( expdir0 + exp + '/' + p, 'rb' ) )
            # save the results
            ds = nc4.Dataset(out_fn, 'a')
            NTx, NPx = ds['lon'].shape
            NT, NP = P['lon'][:,::pstep].shape
            for vn in vlist:
                if vn == 'ot':
                    ds[vn][NTx:] = P[vn][1:]
                elif vn == 'age':
                    ds[vn][NTx:,:] = np.ones((NT-1,1)) * P[vn][::pstep]
                else:
                    ds[vn][NTx:,:] = P[vn][1:,::pstep]
            ds.close()
        counter += 1
    print('  - took %0.1f seconds' % (time.time() - tt0))
    sys.stdout.flush()
