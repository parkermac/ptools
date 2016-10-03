# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 08:56:16 2016

@author: PM5
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

#%% setup input locations

in_dir0 = Ldir['parent'] + 'ptools_output/atlantis/'
in_dir = in_dir0 + 'gridded_polygons/'

out_dir = in_dir0 + 'means_and_fluxes/'
Lfun.make_dir(out_dir, clean=True)

R_in_dir0 = Ldir['parent'] + 'roms/output/salish_2006_4_lp/'

#%% load polygon results

gpoly_dict = pickle.load(open(in_dir + 'gpoly_dict.p', 'rb'))

# and specify z levels
z_dict = {0:0, 1:-5, 2:-25, 3:-50, 4:-100, 5:-150, 6:-350}


#%% find averages and fluxes

dt0 = datetime(2006,1,1)

counter = 0
for ndays in range(3, 363): # 209 is 2006.07.29

    # in /data1/parker/roms/output/salish_2006_4_lp
    # we have f2006.01.04 through 2016.12.29
    # dt0 plus 3 to plus 362 days

    tt0 = time.time()

    dt = dt0 + timedelta(days=ndays)
    f_string = 'f' + dt.strftime('%Y.%m.%d')
    R_in_dir = R_in_dir0 + f_string + '/'
    R_fn = R_in_dir + 'low_passed.nc'
    ds = nc.Dataset(R_fn)

    vv_dict = dict()
    vv_list_full = ['salt', 'temp', 'u', 'v', 'w']
    for vn in vv_list_full:
        vv_dict[vn] = ds[vn][:].squeeze()

    print('Working on day ' + f_string)

    if counter == 0:
        [G, S] = zrfun.get_basic_info(R_fn, getT=False)
        zeta = ds['zeta'][0,:,:]
        z_rho, z_w =  zrfun.get_z(G['h'], zeta, S)
        DA = G['DX'] * G['DY']
        DZ = np.diff(z_w, axis=0)

        DXv = G['DX'][:-1, :] + np.diff(G['DX'], axis=0)/2
        DYu = G['DY'][:, :-1]

        DZu = DZ[:, :, :-1] + np.diff(DZ, axis=2)/2
        DZv = DZ[:, :-1, :] + np.diff(DZ, axis=1)/2

        DAHu = DZu * DYu
        DAHv = DZv * DXv

        z_u = z_rho[:, :, :-1] + np.diff(z_rho, axis=2)/2
        z_v = z_rho[:, :-1, :] + np.diff(z_rho, axis=1)/2

    # intitialize result dicts (key = polygon number)
    trans_dict = dict()
    stvwa_dict = dict()

    for npoly in gpoly_dict.keys():

        print('  npoly = ' + str(npoly))

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
        stvwa_arr = np.zeros((len(z_dict)-1, 5))
        vn_list = ['salt', 'temp']
        for vn in vn_list:
            if print_info:
                print('\n*** ' + vn + ' ***')
            vv = vv_dict[vn] #ds[vn][:].squeeze()
            # 3D array, masked where there is land
            for iz in range(len(z_dict.keys())-1):
                z0 = z_dict[iz + 1] # lower z
                z1 = z_dict[iz] # upper z
                zmask = (z_rho > z0) & (z_rho <= z1)
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
                # check lengths
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
                    print('%d:%d mean = %0.3f, vol = %0.3f km3 (npts=%d)' %
                        (z0, z1, vv_mean, vol/1e9, l1))

                if vn == 'salt':
                    stvwa_arr[iz, 0] = vv_mean
                elif vn == 'temp':
                    stvwa_arr[iz, 1] = vv_mean
                    stvwa_arr[iz, 2] = vol

        # find average w through the bottom of a volume
        vn = 'w'
        if print_info:
            print('\n*** ' + vn + ' ***')
        w = vv_dict['w'] #ds['w'][:].squeeze()
        for iz in range(len(z_dict.keys())-1):
            z0 = z_dict[iz + 1] # layer z
            zmask = (z_rho < z0)
            # 3D bool array, True below z0
            # NOTE: using z_rho here SHOULD mean that our horizontal
            # and vertical boundaries are consistent.
            klev = np.sum(zmask, axis=0)
            # 2D array of vertical indices
            klevp = klev[j_in, i_in]
            # vector of vertical indices
            # same length as j_in and i_in
            #
            this_w = w[klevp, j_in, i_in]
            # vector of w values at the z level, inside the polygon,
            # masked where there is land, but including spurious values
            # where k = 0, same length as klevp, j_in and i_in
            this_wm = np.ma.masked_where(klevp==0, this_w)
            # like this_vv, but now ALSO masked where k = 0
            this_DA = DA[j_in, i_in]
            this_DAm = np.ma.masked_where(this_wm.mask, this_DA)
            # check lengths
            l1 = this_wm.compressed().size
            l2 = this_DAm.compressed().size
            if l1 != l2:
                print('Warning: Result vectors are different lengths')
            if (len(klevp) != len(j_in)) or (len(klevp) != len(i_in)):
                print('Warning: Mask vectors are different lengths')
            # do the integrals
            if l1>0:
                area = this_DAm.sum()
                w_mean = (this_wm * this_DAm).sum() / area
            else:
                area = 0.
                w_mean = 0.
            w_trans = w_mean * area

            stvwa_arr[iz, 3] = w_mean
            stvwa_arr[iz, 4] = area

            if print_info:
                print('%d: w_trans = %0.3f, area = %0.3f km2 (npts=%d)' %
                    (z0, w_trans, area/1e6, l1))

        stvwa_dict[npoly] = stvwa_arr

        # find fluxes through the faces

        u = vv_dict['u']
        v = vv_dict['v']
        # 3D arrays, masked where there is land

        if False: # DEBUGGING
            # test of volume conservation
            wmm = .2
            uda = (u * DYu * DZu).sum(axis=0)
            vda = (v * DXv * DZv).sum(axis=0)
            conv = -(np.diff(uda, axis=1)[1:-1,:] + np.diff(vda, axis=0)[:,1:-1])/DA[1:-1,1:-1]
            w_mean_poly = conv[j_in-1, i_in-1].mean()
            # conv should be the implied vertical velocity at the top of the
            # water column due to lack of mass conservation
            import matplotlib.pyplot as plt
            plt.close('all')
            pth = os.path.abspath('../../LiveOcean/plotting')
            if pth not in sys.path:
                sys.path.append(pth)
            import pfun
            fig = plt.figure(figsize=(18, 10))
            ax = fig.add_subplot(121)
            pfun.add_coast(ax)
            ax.axis([-124, -122, 46.8, 49.2])
            pfun.dar(ax)
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')
            cs = ax.pcolormesh(G['lon_psi'], G['lat_psi'], conv*1000, vmin=-wmm,
                               vmax=wmm)
            fig.colorbar(cs)
            ax.set_title('Convergence w (mm/s)')

            ax = fig.add_subplot(122)
            pfun.add_coast(ax)
            ax.axis([-124, -122, 46.8, 49.2])
            pfun.dar(ax)
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')
            cs = ax.pcolormesh(G['lon_psi'], G['lat_psi'], w[-1,1:-1,1:-1]*1000,
                               vmin=-wmm, vmax=wmm)
            fig.colorbar(cs)
            ax.set_title('Real w (mm/s)')
            plt.show()

        trans_arr = np.zeros((len(z_dict)-1, len(per_dict)))
        for iz in range(len(z_dict.keys())-1):
            z0 = z_dict[iz + 1]
            z1 = z_dict[iz]
            zmask_u = (z_u > z0) & (z_u <= z1)
            zmask_v = (z_v > z0) & (z_v <= z1)
            # 3D bool arrays, true in the layer
            um = np.ma.masked_where(~zmask_u, u)
            vm = np.ma.masked_where(~zmask_v, v)
            # velocity 3D arrays, masked under land,
            # AND outside of the vertical layer
            if print_info:
                print('\n***** TRANSPORT: z range = %d:%d m *****' % (z0, z1))

            for iseg in per_dict.keys():
                per = per_dict[iseg]
                JJ = per[:,0]
                II = per[:,1]
                UV = per[:,2]
                PM = per[:,3]
                # vectors of integers
                JJu = JJ[UV==0]
                IIu = II[UV==0]
                PMu = PM[UV==0]
                JJv = JJ[UV==1]
                IIv = II[UV==1]
                PMv = PM[UV==1]
                # shorter vectors of integers, specific to the u- and v-grids

                this_um = um[:, JJu, IIu]
                this_vm = vm[:, JJv, IIv]

                draft_DAHu = DAHu[:, JJu, IIu]
                draft_DAHv = DAHv[:, JJv, IIv]
                this_DAHu = np.ma.masked_where(this_um.mask, draft_DAHu)
                this_DAHv = np.ma.masked_where(this_vm.mask, draft_DAHv)

                 # check lengths
                l1u = this_um.compressed().size
                l2u = this_DAHu.compressed().size
                if l1u != l2u:
                    print('Warning: U Result vectors are different lengths')

                l1v = this_vm.compressed().size
                l2v = this_DAHv.compressed().size
                if l1v != l2v:
                    print('Warning: V Result vectors are different lengths')

                # do the integrals
                if l1u>0:
                    area_u = this_DAHu.sum()
                    trans_u = (this_um * this_DAHu * PMu).sum()
                else:
                    area_u = 0.
                    trans_u = 0.
                #
                if l1v>0:
                    area_v = this_DAHv.sum()
                    trans_v = (this_vm * this_DAHv * PMv).sum()
                else:
                    area_v = 0.
                    trans_v = 0.

                trans = trans_u + trans_v
                trans_arr[iz,iseg] = trans
                area = area_u + area_v

                # trans_arr is an numpy array with rows = z-levels
                # and columns = transport (m3/s) through polygon faces

                trans_dict[npoly] = trans_arr

                if print_info:
                    print('---- iseg=%d: trans = %0.1f m3/s, area = %0.1f m2 (npts=%d)' %
                        (iseg, trans, area, l1u+l1v))


        if False: # DEBUGGING
            # test volume conservation
            # first pack transports bottom to top
            net_trans = trans_arr.sum(axis=1)
            w_trans_arr = stvwa_arr[:,3] * stvwa_arr[:,4]
            w_trans_arr = w_trans_arr[::-1]
            w_trans_arr = np.concatenate((w_trans_arr, np.zeros(1)))
            net_trans = net_trans[::-1]
            full_net_trans = np.cumsum(net_trans)
            for ii in range(len(w_trans_arr)-1):
                print('w_trans = %0.1f, full_net_trans = %0.1f' %
                    (w_trans_arr[ii+1], full_net_trans[ii]))

    ds.close()

    # save the results for this day

    pickle.dump(stvwa_dict, open(out_dir+f_string+'_stvwa.p', 'wb'))
    pickle.dump(trans_dict, open(out_dir+f_string+'_trans.p', 'wb'))
    counter += 1


    print('  * took %0.1f seconds' % (time.time() - tt0))