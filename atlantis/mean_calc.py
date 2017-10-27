"""
Code to calculate mean values of properties in polygons,
in depth layers.

"""

#%% Imports

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()
import zrfun

import numpy as np
from datetime import datetime, timedelta
import pickle
import netCDF4 as nc
import time

print_info = False
testing = False

#%% setup input locations

whichyear = 2005
if Ldir['env'] == 'pm_mac': # mac version
    if whichyear == 2006:
        R_in_dir0 = Ldir['parent'] + 'roms/output/salish_2006_4_lp/'
        out_dir0 = Ldir['parent'] + 'ptools_output/atlantis_mac_2006/'
    elif whichyear == 2005:
        R_in_dir0 = Ldir['parent'] + 'roms/output/salish_2005_1_lp/'
        out_dir0 = Ldir['parent'] + 'ptools_output/atlantis_mac_2005/'
elif Ldir['env'] == 'fjord': # fjord version
    if whichyear == 2006:
        R_in_dir0 = '/boildat1/parker/roms/output/salish_2006_4_lp/'
        out_dir0 = Ldir['parent'] + 'ptools_output/atlantis_fjord_2006/'
    elif whichyear == 2005:
        R_in_dir0 = '/boildat1/parker/roms/output/salish_2005_1_lp/'
        out_dir0 = Ldir['parent'] + 'ptools_output/atlantis_fjord_2005/'

in_dir = out_dir0 + 'gridded_polygons/'

out_dir = out_dir0 + 'means/'
Lfun.make_dir(out_dir, clean=True)

# load polygon results

gpoly_dict = pickle.load(open(in_dir + 'gpoly_dict.p', 'rb'))

# and specify z levels
z_dict = {0:5, 1:-5, 2:-25, 3:-50, 4:-100, 5:-150, 6:-350}
NLAY = len(z_dict) - 1

#%% find averages

dt0 = datetime(whichyear,1,1)

if Ldir['env'] == 'pm_mac':
    # have 2006.07.01-31 = days 181 to 211
    # big convergence errors for 7/29, 7/30 = 209, 210
    day_list = [208, 209] #range(181,211+1)
elif Ldir['env'] == 'fjord':
    # in /data1/parker/roms/output/salish_2006_4_lp
    # we have f2006.01.04 through 2016.12.29
    # = days 3 to 362
    day_list = range(3, 363)

counter = 0
for nday in day_list:

    tt0 = time.time()

    dt = dt0 + timedelta(days=nday)
    f_string = 'f' + dt.strftime('%Y.%m.%d')
    R_in_dir = R_in_dir0 + f_string + '/'
    R_fn = R_in_dir + 'low_passed.nc'
        
    ds = nc.Dataset(R_fn)

    vv_dict = dict()
    vv_list_full = ['salt', 'temp']
    for vn in vv_list_full:
        vv_dict[vn] = ds[vn][:].squeeze()

    print('\nWorking on day ' + f_string)

    if counter == 0:
        G = zrfun.get_basic_info(R_fn, only_G=True)
        S = zrfun.get_basic_info(R_fn, only_S=True)
        zeta = ds['zeta'][0,:,:]
        z_rho, z_w =  zrfun.get_z(G['h'], zeta, S)
        DA = G['DX'] * G['DY']
        DZ = np.diff(z_w, axis=0)

    # intitialize result dicts (key = polygon number)
    stv_dict = dict()

    if testing == True:
        npoly_list = range(2)
    else:
        npoly_list = gpoly_dict.keys()
        
    for npoly in npoly_list:

        if print_info:
            print('\n>>>> npoly = ' + str(npoly) + ' <<<<')

        # we have two objects associated with a given polygon,
        #  * per_dict[iseg] has arrays of boundary information, and
        #  * ji_rho_in is an array of indices of interior points
        per_dict = gpoly_dict[npoly]['per_dict']
        ji_rho_in = gpoly_dict[npoly]['ji_rho_in']

        pmask = np.ones_like(DA) == 0
        # start with pmask False everywhere
        j_in = ji_rho_in[:,0]
        i_in = ji_rho_in[:,1]
        pmask[j_in, i_in] = True
        # pmask is True inside the polygon

        DAm = np.ma.masked_where(~pmask, DA)
        # 2D array of horizontal areas,
        # masked outside of the polygon

        DVm = DZ * DAm
        # 3D array of cell volumes,
        # masked outside of the polygon

        # find average values inside the volume
        stv_arr = np.zeros((NLAY, 3))
        vn_list = ['salt', 'temp']
        for vn in vn_list:
            if print_info:
                print(' *** ' + vn + ' ***')
            vv = vv_dict[vn]
            # 3D array, masked where there is land
            for iz in range(NLAY):
                z_lower = z_dict[iz + 1] # lower z
                z_upper = z_dict[iz] # upper z
                zmask = (z_rho > z_lower) & (z_rho <= z_upper)
                # 3D bool array, true in the layer
                #
                vvDVm = vv*DVm
                # 3D array
                # masked where there is land AND outside of the polygon
                #
                this_vvDVm = vvDVm[zmask]
                draft_DVm = DVm[zmask]
                # ensure that the mask on the volume is identical to that
                # on the data
                this_DVm = np.ma.masked_where(this_vvDVm.mask, draft_DVm)
                # numpy vectors, all with the same mask
                #
                # check lengths [.compressed().size counts only the non-masked values]
                l1 = this_DVm.compressed().size
                l2 = this_vvDVm.compressed().size
                if l1 != l2:
                    print('Warning: Result vectors are different lengths')
                # do the integrals
                if l1>0:
                    vol = this_DVm.sum()
                    vv_mean = this_vvDVm.sum() / vol
                else:
                    vol = 0.
                    vv_mean = np.nan
                    
                if print_info:
                    print(' %4d:%4d mean = %6.2f, vol = %6.2f km3 (npts=%6d)' %
                        (z_lower, z_upper, vv_mean, vol/1e9, l1))

                if vn == 'salt':
                    stv_arr[iz, 0] = vv_mean
                elif vn == 'temp':
                    stv_arr[iz, 1] = vv_mean
                    stv_arr[iz, 2] = vol

        stv_dict[npoly] = stv_arr

    ds.close()

    # save the results for this day

    pickle.dump(stv_dict, open(out_dir+f_string+'_stv.p', 'wb'))
    counter += 1

    print('  Took %0.1f seconds' % (time.time() - tt0))