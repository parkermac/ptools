"""
Code to calculate mean fluxes of properties through polygon
faces, in depth layers.

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
import matplotlib.pyplot as plt

testing = False

if testing == True:
    plt.close('all')

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

out_dir = out_dir0 + 'fluxes/'
Lfun.make_dir(out_dir, clean=True)

# load polygon results

gpoly_dict = pickle.load(open(in_dir + 'gpoly_dict.p', 'rb'))
shared_faces_dict = pickle.load(open(in_dir + 'shared_faces.p', 'rb'))

# and specify z levels
z_dict = {0:5, 1:-5, 2:-25, 3:-50, 4:-100, 5:-150, 6:-350}
NLAY = len(z_dict) - 1

#%% find fluxes

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
    
    # specify ROMS file to work on
    dt0 = datetime(2006,1,1)
    dt = dt0 + timedelta(days=nday)
    f_string = 'f' + dt.strftime('%Y.%m.%d')
    print('\nWorking on day %s (nday = %3d)' % (f_string, nday))
    R_in_dir = R_in_dir0 + f_string + '/'
    R_fn = R_in_dir + 'low_passed.nc'
    ds = nc.Dataset(R_fn)

    u = ds['u'][:].squeeze()
    v = ds['v'][:].squeeze()
    w0 = ds['w'][0,-1,:,:].squeeze()

    G = zrfun.get_basic_info(R_fn, only_G=True)
    S = zrfun.get_basic_info(R_fn, only_S=True)
    zeta = ds['zeta'][0,:,:]
    z_rho, z_w =  zrfun.get_z(G['h'], zeta, S)
    DA = G['DX'] * G['DY']
    DAm = np.ma.masked_where(zeta.mask, DA)
    DZ = np.diff(z_w, axis=0)
    
    # make versions of DA masked for the z layers
    zDAm = DAm * np.ones((NLAY, 1, 1))
    for nlay in range(NLAY):
        z_lower = z_dict[nlay + 1] # lower z
        zmask = (-G['h'] >= z_lower)
        draft_DA = DA.copy()
        draft_DA[zmask] = 0.
        this_DAm = np.ma.masked_where(zeta.mask, draft_DA)
        zDAm[nlay,:,:] = this_DAm
            
    # Z on u and v grids
    Zu = z_rho[:, :, :-1] + np.diff(z_rho, axis=2)/2
    Zv = z_rho[:, :-1, :] + np.diff(z_rho, axis=1)/2
    zmu_dict = dict()
    zmv_dict = dict()
    for nlay in range(NLAY):
        z_lower = z_dict[nlay + 1] # lower z
        z_upper = z_dict[nlay] # upper z
        zmu_dict[nlay] = (Zu > z_lower) & (Zu <= z_upper)            
        zmv_dict[nlay] = (Zv > z_lower) & (Zv <= z_upper)
    layu_dict = dict()
    layv_dict = dict()
    
    # DZ on u and v grids
    DZu = DZ[:, :, :-1] + np.diff(DZ, axis=2)/2
    DZv = DZ[:, :-1, :] + np.diff(DZ, axis=1)/2
    # DX and DY on u and v grids
    DYu = G['DY'][:, :-1]
    DXv = G['DX'][:-1, :] + np.diff(G['DX'], axis=0)/2
    # cell areas for u and v grid box faces
    DAHu = DZu * DYu
    DAHv = DZv * DXv
    
    # Initialize arrays to store transport data.
    face_trans_dict = dict()
    poly_conv_dict = dict()

    face_ztrans_dict = dict()
    poly_zconv_dict = dict()
    
    # calculate convergence in each polygon
    counter = 0
    NPOLY = len(gpoly_dict)
    for npoly in range(NPOLY):
    
        #print('  npoly = ' + str(npoly))

        # we have two objects associated with a given polygon,
        #  * per_dict[iseg] has arrays of boundary information, and
        #  * ji_rho_in is an array of indices of interior points
        per_dict = gpoly_dict[npoly]['per_dict']
        ji_rho_in = gpoly_dict[npoly]['ji_rho_in']
    
        j_in = ji_rho_in[:,0]
        i_in = ji_rho_in[:,1]
    
        # find fluxes through the faces
        NFACE = len(per_dict)                
        face_trans_arr = np.zeros(NFACE)
        face_area_arr = np.zeros(NFACE)
        
        face_ztrans_arr = np.zeros((NFACE, NLAY))
        face_zarea_arr = np.zeros((NFACE, NLAY))

        for nface in range(NFACE):
            per = per_dict[nface]                
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
            PMu = PMu.reshape((1, len(PMu)))
            PMv = PMv.reshape((1, len(PMv)))
            this_u = u[:, JJu, IIu]
            this_v = v[:, JJv, IIv]
                        
            draft_DAHu = DAHu[:, JJu, IIu]
            draft_DAHv = DAHv[:, JJv, IIv]
            this_DAHu = np.ma.masked_where(this_u.mask, draft_DAHu)
            this_DAHv = np.ma.masked_where(this_v.mask, draft_DAHv)
             # check lengths
            l1u = this_u.compressed().size
            l2u = this_DAHu.compressed().size
            if l1u != l2u:
                print('Warning: U Result vectors are different lengths')
            l1v = this_v.compressed().size
            l2v = this_DAHv.compressed().size
            if l1v != l2v:
                print('Warning: V Result vectors are different lengths')
            # do the integrals
            if l1u>0:
                area_u = this_DAHu.sum()
                trans_u = (this_u * this_DAHu * PMu).sum()
            else:
                area_u = 0.
                trans_u = 0.
            if l1v>0:
                area_v = this_DAHv.sum()
                trans_v = (this_v * this_DAHv * PMv).sum()
            else:
                area_v = 0.
                trans_v = 0.

            face_trans = trans_u + trans_v
            face_trans_arr[nface] = face_trans
            face_area = area_u + area_v
            face_area_arr[nface] = face_area
        
            # store results for later
            face_trans_dict[(npoly, nface)] = (face_trans, face_area)
            
            # START z level code #############################################
            
            # now do the same thing but divvying up into z levels
            for nlay in range(NLAY):
                this_zmu = zmu_dict[nlay][:, JJu, IIu]
                this_zmv = zmv_dict[nlay][:, JJv, IIv]
                
                this_zu = this_u[this_zmu]
                this_zv = this_v[this_zmv]
                                
                draft_zDAHu = this_DAHu[this_zmu]
                draft_zDAHv = this_DAHv[this_zmv]
                
                this_zDAHu = np.ma.masked_where(this_zu.mask, draft_zDAHu)
                this_zDAHv = np.ma.masked_where(this_zv.mask, draft_zDAHv)
            
                 # check lengths
                l1zu = this_zu.compressed().size
                l2zu = this_zDAHu.compressed().size
                if l1zu != l2zu:
                    print('Warning: ZU Result vectors are different lengths')
                l1zv = this_zv.compressed().size
                l2zv = this_zDAHv.compressed().size
                if l1zv != l2zv:
                    print('Warning: ZV Result vectors are different lengths')
                # do the integrals
                if l1zu>0:
                    area_zu = this_zDAHu.sum()
                    trans_zu = (this_zu * this_zDAHu * PMu[0,0]).sum()
                else:
                    area_zu = 0.
                    trans_zu = 0.
                if l1zv>0:
                    area_zv = this_zDAHv.sum()
                    trans_zv = (this_zv * this_zDAHv * PMv[0,0]).sum()
                else:
                    area_zv = 0.
                    trans_zv = 0.

                face_ztrans = trans_zu + trans_zv
                face_ztrans_arr[nface, nlay] = face_ztrans
                face_zarea = area_zu + area_zv
                face_zarea_arr[nface, nlay] = face_zarea
                # Do we use the _arr fields?  I think not, but it doesn't matter.            
        
                # store results for later
                face_ztrans_dict[(npoly, nface, nlay)] = (face_ztrans, face_zarea)
            
            # END z level code ###############################################
            
        # check results of z level code
        # RESULT it works perfectly (to within roundoff error I think)
        fzt = 0.
        fza = 0.
        for nlay in range(NLAY):
            fzt += face_ztrans_dict[(npoly, nface, nlay)][0]
            fza += face_ztrans_dict[(npoly, nface, nlay)][1]
        if np.abs(fzt - face_trans_dict[(npoly, nface)][0]) > .001:
            print('npoly=%d nface=%d transport error' % (npoly, nface))
            print('fzt=%0.5f ft=%0.5f' % (fzt, face_trans_dict[(npoly, nface)][0]))
        if np.abs(fza - face_trans_dict[(npoly, nface)][1]) > .001:
            print('npoly=%d nface=%d area error' % (npoly, nface))
            print('fza=%0.5f fa=%0.5f' % (fza, face_trans_dict[(npoly, nface)][1]))

        poly_zarea = np.zeros(NLAY)
        for nlay in range(NLAY):
            this_zarea = zDAm[nlay,j_in,i_in].sum()
            try:
                if this_zarea.mask == True:
                    this_zarea = 0.
            except AttributeError:
                pass            
            poly_zarea[nlay] = this_zarea
            
        poly_area = DAm[j_in,i_in].sum()
        try:
            if poly_area.mask == True:
                poly_area = 0.
        except AttributeError:
            pass
        net_conv = face_trans_arr.sum()       
        if (poly_area > 0):
            poly_mean_w = net_conv/poly_area
        else:
            poly_mean_w = 0.0
        
        # store results for later
        net_face_area = face_area_arr.sum() 
        poly_conv_dict[npoly] = (net_conv, poly_area, poly_zarea, net_face_area, NFACE)
    
        counter += 1

    ds.close()
    
    # save originals
    orig_poly_conv_dict = poly_conv_dict.copy()
    orig_face_trans_dict = face_trans_dict.copy()
    orig_face_ztrans_dict = face_ztrans_dict.copy()

    # Next try to adjust all polygons to have conv = 0
    NITER = 400
    for iii in range(NITER):
    
        new_poly_conv_dict = poly_conv_dict.copy()
        for npoly in range(NPOLY):
            net_conv, poly_area, poly_zarea, net_face_area, NFACE = new_poly_conv_dict[npoly]
            new_poly_conv_dict[npoly] = (0.0, poly_area, poly_zarea, net_face_area, NFACE)
    
        new_face_trans_dict = face_trans_dict.copy()

        for npoly in range(NPOLY):
            net_conv, poly_area, poly_zarea, net_face_area, NFACE = poly_conv_dict[npoly]
            new_net_conv, poly_area, poly_zarea, net_face_area, NFACE = new_poly_conv_dict[npoly]
            dconv = new_net_conv - net_conv
    
            if net_face_area != 0.0:
                dconv_a = dconv / net_face_area
            else:
                dconv_a = 0.0
    
            for nface in range(NFACE):
                face_trans, face_area = face_trans_dict[(npoly, nface)]
                if face_trans != 0.0 and face_area != 0.0:
                    new_face_trans_dict[(npoly, nface)] = (face_trans + dconv_a*face_area, face_area)
                else:
                    pass # keep original values

        new_face_trans_dict_copy = new_face_trans_dict.copy()
        for npoly in range(NPOLY):
            net_conv, poly_area, poly_zarea, net_face_area, NFACE = poly_conv_dict[npoly]
            for nface in range(NFACE):
                try:
                    new_face_trans, face_area = new_face_trans_dict[(npoly, nface)]
                    ipoly, iface = shared_faces_dict[(npoly, nface)]
                    new_facing_trans, facing_area = new_face_trans_dict_copy[(ipoly, iface)]
                    fact = (new_face_trans + new_facing_trans)/2
                    new_face_trans_dict[(npoly, nface)] = (new_face_trans - fact, face_area)
                except KeyError:
                    # presumably this face does not have a match on another polygon
                    pass

        for npoly in range(NPOLY):
            net_conv, poly_area, poly_zarea, net_face_area, NFACE = poly_conv_dict[npoly]
            shelf = []
            for nface in range(NFACE):
                new_face_trans, face_area = new_face_trans_dict[(npoly, nface)]
                shelf.append(new_face_trans)
            new_conv = np.array(shelf).sum()
            new_poly_conv_dict[npoly] = (new_conv, poly_area, poly_zarea, net_face_area, NFACE)
            
        face_trans_dict = new_face_trans_dict.copy()
        poly_conv_dict = new_poly_conv_dict.copy()
    
    # finally add the adjustments to the transports in z levels:
    face_ztrans_dict = dict()
    for npoly in range(NPOLY):
        net_conv, poly_area, poly_zarea, net_face_area, NFACE = poly_conv_dict[npoly]
        for nface in range(NFACE):
            face_trans, face_area = face_trans_dict[(npoly, nface)]
            orig_face_trans, face_area = orig_face_trans_dict[(npoly, nface)]
            DQ = face_trans - orig_face_trans
            for nlay in range(NLAY):
                (face_ztrans, face_zarea) = orig_face_ztrans_dict[(npoly, nface, nlay)]
                if face_area > 0:
                    adj = DQ * face_zarea / face_area
                else:
                    adj = 0.
                face_ztrans_dict[(npoly, nface, nlay)] = (face_ztrans + adj, face_zarea)

    # calculate w at interfaces and check that they are reasonable
    # RESULT: they seem good
    wz_dict = dict()
    # each entry in wz_dict is a tuple of (w, area) at the TOP of a layer
    print_info = False
    for npoly in range(NPOLY):
        if print_info:
            print('\nnpoly = %d' % (npoly))
        cz = np.zeros(NLAY)
        net_conv, poly_area, poly_zarea, net_face_area, NFACE = poly_conv_dict[npoly]
        for nface in range(NFACE):
            face_trans, face_area = face_trans_dict[(npoly, nface)]
            for nlay in range(NLAY):
                (face_ztrans, face_zarea) = face_ztrans_dict[(npoly, nface, nlay)]
                cz[nlay] += face_ztrans
        czr = cz[::-1] # packed bottom to top
        cczr = np.cumsum(czr)
        ccz = cczr[::-1]
        # ccz is the vertical transport through the UPPER boundary of each layer
        # packed top to bottom like all other z-layer variables
        for nlay in range(NLAY):
            if nlay == 0: # get the associated horizontal area
                this_zarea = poly_area
            else:
                this_zarea = poly_zarea[nlay-1]
            try: # accout for the few cases where an area might be a masked constant
                if this_zarea.mask == True:
                    this_zarea = 0.
            except AttributeError:
                pass
            # calculate the vertical velocity through the layer
            # note that the transport throught the deepest level (-350 m) is ZERO
            # and the uppermost w is also ~zero because of our iterative correction
            if (this_zarea == 0.):
                wz_dict[(npoly, nlay)] = (0., 0.)
            else:
                wz_dict[(npoly, nlay)] = (ccz[nlay] / this_zarea, this_zarea)
            if print_info:
                print(' z = %4d w = %10.1f (mm/hour)' %
                    (z_dict[nlay], 3600*1000*wz_dict[(npoly, nlay)][0]))

    # save the results for this day
    pickle.dump(face_trans_dict, open(out_dir+f_string+'_face_trans.p', 'wb'))
    pickle.dump(face_ztrans_dict, open(out_dir+f_string+'_face_ztrans.p', 'wb'))
    pickle.dump(poly_conv_dict, open(out_dir+f_string+'_poly_conv.p', 'wb'))
    pickle.dump(wz_dict, open(out_dir+f_string+'_wz.p', 'wb'))

    # check that the net convergence is still small
    # RESULT: it is very,very small
    zconv = []
    pconv = []
    for npoly in range(NPOLY):
        cc = 0.
        net_conv, poly_area, poly_zarea, net_face_area, NFACE = poly_conv_dict[npoly]
        for nface in range(NFACE):
            face_trans, face_area = face_trans_dict[(npoly, nface)]
            for nlay in range(NLAY):
                (face_ztrans, face_zarea) = face_ztrans_dict[(npoly, nface, nlay)]
                cc += face_ztrans
        zconv.append(cc)
        pconv.append(net_conv)
    print('  max convergence error = %g (m3/s)' % (np.abs(np.array(pconv) - np.array(zconv)).max()))
                                
    # check size of w at the free surface
    w_arr = np.zeros(NPOLY)
    for k in range(NPOLY):
        conv, poly_area, poly_zarea, net_face_area, NFACE = poly_conv_dict[k]
        if poly_area > 0:
            w_arr[k] = conv / poly_area
        else:
            w_arr[k] = 0.
    print('  max w = %0.5f (mm/hour)' % (3600 * 1000 * np.abs(w_arr.max())))
    
    # report on calculation time
    print('  Took %0.1f seconds' % (time.time() - tt0))

    # a plot of the net convergence in all polygons, before and after
    # the iterative correction
    if (testing == True):
        new_conv = []
        orig_conv = []
        ipoly = []
        for k in poly_conv_dict.keys():
            new_conv.append(poly_conv_dict[k][0])
            orig_conv.append(orig_poly_conv_dict[k][0])
            ipoly.append(k)
            
        fig = plt.figure(figsize=(10,6))
        ax = fig.add_subplot(111)
        ax.plot(ipoly, orig_conv, '*r', ipoly, new_conv, 'ob')
        ax.set_title(f_string)
        ax.set_xlabel('Polygon Number')
        ax.set_ylabel('Convergence (m3/s)')
        plt.show()
                