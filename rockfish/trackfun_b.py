"""
Functions for particle tracking.
"""
# setup
import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:sys.path.append(alp)
import zfun
import zrfun
import numpy as np
import netCDF4 as nc4
from datetime import datetime, timedelta

def get_tracks(fn_list, plon0, plat0, pcs0, dir_tag, method,
                surface, turb, ndiv, windage, bound='stop', dep_range=None):
    '''
    Create tracks of particles starting from the locations at (plon0, plat0, pcs0)
    
    Inputs:
    fn_list - array of dataset file locations
    plon0/plat0/pcs0 - arrays of position (longitude/latitude/proportion of total depth)
    dir_tag - string determining direction of tracking ('forward' or 'reverse')
    method - string determining interpolation level ('rk2' or 'rk4')
    surface - Boolean, True traps particles to the surface
    turb - Boolean, True adds random vertical walk
    ndiv - int determining number of divisions to make between saves for the integration
    windage - float (0 to 1) determining proportion of velocity to add from surface wind
    bound - string determing boundary velocity conditions 
        'stop' (default) prevents particles from moving when a gridcell contains land
        'reflect' causes particles to move away from gridcells that contain land
    dep_range - tuple containing arrays (minimum depth, maximum depth) for all particles
    
    Outputs:
    P - dictionary containing variables in plist_main
        each variable has rows of time and columns of particles
        except the array 'ot', seconds since 1970
    G - dictionary containing grid and bathymetry information
    S - dictionary containing water surface information
    '''

# Bartos - removed on land particle removal, so set these to plon, not plonA
    plon = plon0.copy()
    plat = plat0.copy()
    pcs = pcs0.copy()

    # get basic info
    G = zrfun.get_basic_info(fn_list[0], only_G=True)
    S = zrfun.get_basic_info(fn_list[0], only_S=True)

    # get time vector of history files
    NT = len(fn_list)
    rot = np.nan * np.ones(NT)
    counter = 0
    for fn in fn_list:
        ds = nc4.Dataset(fn)
        rot[counter] = ds.variables['ocean_time'][:].squeeze()
        counter += 1
        ds.close

    delta_t = rot[1] - rot[0] # seconds between saves

    # this is how we track backwards in time
    if dir_tag == 'reverse':
        delta_t = -delta_t
        fn_list = fn_list[::-1]

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

    # these lists are used internally to get other variables as needed
    vn_list_vel = ['u','v','w']
    vn_list_zh = ['zeta','h']
    vn_list_wind = ['Uwind','Vwind']
    vn_list_other = ['salt', 'temp', 'zeta', 'h', 'u', 'v', 'w']

    # Step through times.
# Bartos - removed on land particle removal
    counter = 0
    nrot = len(rot)
    for pot in rot[:-1]:

        if np.mod(counter,24) == 0:
            print(' - time %d out of %d' % (counter, nrot))
            sys.stdout.flush()

        # get time indices
        it0, it1, frt = zfun.get_interpolant(
                np.array(pot), rot, extrap_nan=False)

        # get the velocity zeta, and h at all points
        ds0 = nc4.Dataset(fn_list[it0[0]])
        ds1 = nc4.Dataset(fn_list[it1[0]])

        if counter == 0:

            # create result arrays
            NP = len(plon)
            P = dict()
            # plist main is what ends up written to output
            plist_main = ['lon', 'lat', 'cs', 'ot', 'z', 'zeta', 'zbot',
                          'salt', 'temp', 'u', 'v', 'w', 'Uwind', 'Vwind', 'h']
            for vn in plist_main:
                P[vn] = np.nan * np.ones((NT,NP))

            # write positions to the results arrays
            P['lon'][it0,:] = plon
            P['lat'][it0,:] = plat
            if surface == True:
                pcs[:] = S['Cs_r'][-1]
            P['cs'][it0,:] = pcs

            P = get_properties(vn_list_other, ds0, it0, P, plon, plat, pcs, R, surface)

        delt = delta_t/ndiv

        for nd in range(ndiv):

            fr0 = nd/ndiv
            fr1 = (nd + 1)/ndiv
            frmid = (fr0 + fr1)/2

            if method == 'rk4':
                # RK4 integration

                V0, ZH0 = get_vel(vn_list_vel, vn_list_zh, ds0, ds1, 
                                  plon, plat, pcs, R, fr0, surface, bound=bound)

                plon1, plat1, pcs1 = update_position(V0, ZH0, S, delt/2,
                                                     plon, plat, pcs, surface)
                V1, ZH1 = get_vel(vn_list_vel, vn_list_zh, ds0, ds1, 
                                  plon1, plat1, pcs1, R, frmid, surface, bound=bound)

                plon2, plat2, pcs2 = update_position(V1, ZH1, S, delt/2,
                                                     plon, plat, pcs, surface)
                V2, ZH2 = get_vel(vn_list_vel, vn_list_zh, ds0, ds1, 
                                  plon2, plat2, pcs2, R, frmid, surface, bound=bound)

                plon3, plat3, pcs3 = update_position(V2, ZH2, S, delt,
                                                     plon, plat, pcs, surface)
                V3, ZH3 = get_vel(vn_list_vel, vn_list_zh, ds0, ds1, 
                                  plon3, plat3, pcs3, R, fr1, surface, bound=bound)

                # add windage, calculated from the middle time
                if (surface == True) and (windage > 0):
                    Vwind = get_wind(vn_list_wind, ds0, ds1, plon, plat, pcs, R, frmid, surface)
                    Vwind3 = np.concatenate((windage*Vwind, np.zeros((NP,1))), axis=1)
                else:
                    Vwind3 = np.zeros((NP,3))
                
                plon, plat, pcs = update_position((V0 + 2*V1 + 2*V2 + V3)/6 + Vwind3,
                                                  (ZH0 + 2*ZH1 + 2*ZH2 + ZH3)/6,
                                                  S, delt, plon, plat, pcs, surface)
                                                  
# Bartos - begin turbulence edit
                
                # add turbulence in two distinct timesteps
                if turb == True:
                    # pull values of VdAKs and add up to 3-dimensions
                    VdAKs = get_dAKs(vn_list_zh, ds0, ds1, plon, plat, pcs, R, S, frmid, surface)
                    VdAKs3 = np.concatenate((np.zeros((NP,2)), VdAKs[:,np.newaxis]), axis=1)
                    
                    # update position with 1/2 of AKs gradient
                    plon, plat, pcs = update_position(VdAKs3/2, (ZH0 + 2*ZH1 + 2*ZH2 + ZH3)/6,
                                                  S, delt, plon, plat, pcs, surface)
                    
                    # update position with rest of turbulence
                    Vturb = get_turb(ds0, ds1, VdAKs, delta_t, plon, plat, pcs, R, frmid, surface)
                    Vturb3 = np.concatenate((np.zeros((NP,2)), Vturb[:,np.newaxis]), axis=1)
                    
                    plon, plat, pcs = update_position(Vturb3, (ZH0 + 2*ZH1 + 2*ZH2 + ZH3)/6,
                                                  S, delt, plon, plat, pcs, surface)
                
 # Bartos - end edit
                
            elif method == 'rk2':
                # RK2 integration
                V0, ZH0 = get_vel(vn_list_vel, vn_list_zh, ds0, ds1,
                                  plon, plat, pcs, R, fr0, surface, bound=bound)

                plon1, plat1, pcs1 = update_position(V0, ZH0, S, delt/2,
                                                     plon, plat, pcs, surface)
                V1, ZH1 = get_vel(vn_list_vel, vn_list_zh, ds0, ds1,
                                  plon1, plat1, pcs1, R, frmid, surface, bound=bound)

                # add windage, calculated from the middle time
                if (surface == True) and (windage > 0):
                    Vwind = get_wind(vn_list_wind, ds0, ds1, plon, plat, pcs, R, frmid, surface)
                    Vwind3 = np.concatenate((windage*Vwind,np.zeros((NP,1))),axis=1)
                else:
                    Vwind3 = np.zeros((NP,3))

                plon, plat, pcs = update_position(V1 + Vwind3, ZH1, S, delt,
                                                  plon, plat, pcs, surface)
                                                  
# Bartos - begin turbulence edit
                
                # add turbulence in two distinct timesteps
                if turb == True:
                    # pull values of VdAKs and add up to 3-dimensions
                    VdAKs = get_dAKs(vn_list_zh, ds0, ds1, plon, plat, pcs, R, S, frmid, surface)
                    VdAKs3 = np.concatenate((np.zeros((NP,2)), VdAKs[:,np.newaxis]), axis=1)
                    
                    # update position with 1/2 of AKs gradient
                    plon, plat, pcs = update_position(VdAKs3/2, ZH1, S, delt, plon, 
                                                  plat, pcs, surface)
                    
                    # update position with rest of turbulence
                    Vturb = get_turb(ds0, ds1, VdAKs, delta_t, plon, plat, pcs, R, frmid, surface)
                    Vturb3 = np.concatenate((np.zeros((NP,2)), Vturb[:,np.newaxis]), axis=1)
                    
                    plon, plat, pcs = update_position(Vturb3, ZH1, S, delt, plon, 
                                                  plat, pcs, surface)
                                                  
        # write positions to the results arrays
        P['lon'][it1,:] = plon
        P['lat'][it1,:] = plat
        if surface == True:
            pcs[:] = S['Cs_r'][-1]
        P['cs'][it1,:] = pcs
        P = get_properties(vn_list_other, ds1, it1, P, plon, plat, pcs, R, surface)

# Bartos - begin maximum and minimum depth edit

        if dep_range != None:
            dep_max = dep_range[1]
            # mask of depths below dep_max
            dep_mask = P['z'][it1,:]<dep_max
            dep_mask = dep_mask.squeeze()
            # total depth
            dep_tot = P['zeta'][it1,dep_mask] + P['h'][it1,dep_mask]
            # replace masked depths with dep_max
            P['z'][it1,dep_mask] = dep_max[dep_mask]
            P['cs'][it1,dep_mask] = dep_max[dep_mask] / dep_tot
            pcs[dep_mask] = dep_max[dep_mask] / dep_tot
            
            dep_min = dep_range[0]
            # mask of depths above dep_min
            dep_mask = P['z'][it1,:]>dep_min
            dep_mask = dep_mask.squeeze()
            # total depth
            dep_tot = P['zeta'][it1,dep_mask] + P['h'][it1,dep_mask]
            # replace masked depths with dep_min
            P['z'][it1,dep_mask] = dep_min[dep_mask]
            P['cs'][it1,dep_mask] = dep_min[dep_mask] / dep_tot
            pcs[dep_mask] = dep_max[dep_mask] / dep_tot

# Bartos - end edits

        ds0.close()
        ds1.close()

        counter += 1

    # by doing this the points are going forward in time
    if dir_tag == 'reverse':
        for vn in plist_main:
            P[vn] = P[vn][::-1,:]

    # and save the time vector (seconds in whatever the model reports)
    P['ot'] = rot

    return P, G, S

def update_position(V, ZH, S, delta_t, plon, plat, pcs, surface):
    # find the new position
    Plon = plon.copy()
    Plat = plat.copy()
    Pcs = pcs.copy()
    dX_m = V*delta_t
    per_m = zfun.earth_rad(Plat)
    clat = np.cos(np.pi*Plat/180.)
    pdx_deg = (180./np.pi)*dX_m[:,0]/(per_m*clat)
    pdy_deg = (180./np.pi)*dX_m[:,1]/per_m
    H = ZH.sum(axis=1)
    pdz_s = dX_m[:,2]/H
    Plon += pdx_deg
    Plat += pdy_deg
    if surface == False:
        Pcs_orig = Pcs.copy()
        Pcs += pdz_s
        # enforce limits on cs
        mask = np.isnan(Pcs)
        Pcs[mask] = Pcs_orig[mask]
        Pcs[Pcs < S['Cs_r'][0]] = S['Cs_r'][0]
        Pcs[Pcs > S['Cs_r'][-1]] = S['Cs_r'][-1]
    else:
        Pcs[:] = S['Cs_r'][-1]

    return Plon, Plat, Pcs

def get_vel(vn_list_vel, vn_list_zh, ds0, ds1, plon, plat, pcs, R, frac, surface, bound='stop'):
    # get the velocity, zeta, and h at all points, at an arbitrary
    # time between two saves
    # "frac" is the fraction of the way between the times of ds0 and ds1
    # 0 <= frac <= 1
    V0 = get_V(vn_list_vel, ds0, plon, plat, pcs, R, surface, bound=bound)
    V1 = get_V(vn_list_vel, ds1, plon, plat, pcs, R, surface, bound=bound)
    V0[np.isnan(V0)] = 0.0
    V1[np.isnan(V1)] = 0.0
    ZH0 = get_V(vn_list_zh, ds0, plon, plat, pcs, R, surface)
    ZH1 = get_V(vn_list_zh, ds1, plon, plat, pcs, R, surface)
    V = (1 - frac)*V0 + frac*V1
    ZH = (1 - frac)*ZH0 + frac*ZH1

    return V, ZH

def get_wind(vn_list_wind, ds0, ds1, plon, plat, pcs, R, frac, surface):
    # get the wind velocity at an arbitrary
    # time between two saves
    # "frac" is the fraction of the way between the times of ds0 and ds1
    # 0 <= frac <= 1
    V0 = get_V(vn_list_wind, ds0, plon, plat, pcs, R, surface)
    V1 = get_V(vn_list_wind, ds1, plon, plat, pcs, R, surface)
    V0[np.isnan(V0)] = 0.0
    V1[np.isnan(V1)] = 0.0
    V = (1 - frac)*V0 + frac*V1

    return V

# Bartos - begin edit to add functions: 
#          get_dAKS and get_turb for turbulence, change position to 
#          move particles near land

def get_dAKs(vn_list_zh, ds0, ds1, plon, plat, pcs, R, S, frac, surface):
    # create diffusivity gradient for turbulence calculation
    
    # first time step
    ZH0 = get_V(vn_list_zh, ds0, plon, plat, pcs, R, surface)
    ZH0[np.isnan(ZH0)] = 0
    dpcs0 = 1/(ZH0[:,0] + ZH0[:,1]) # change in pcs for a total of a 2m difference
    
    #     upper variables
    pcs0u = pcs-dpcs0
    pcs0u[pcs0u > S['Cs_r'][-1]] = S['Cs_r'][-1]
    AKs0u = get_V(['AKs',], ds0, plon, plat, pcs0u, R, surface) # diffusivity
    z0u = (pcs0u)*(ZH0[:,0]+ZH0[:,1]) # depth = pcs * full-depth
    
    #     lower variables
    pcs0b = pcs+dpcs0
    pcs0b[pcs0b < S['Cs_r'][0]] = S['Cs_r'][0]
    AKs0b = get_V(['AKs',], ds0, plon, plat, pcs0b, R, surface) # diffusivity
    z0b = (pcs0b)*(ZH0[:,0]+ZH0[:,1]) # depth = pcs * full-depth
    
    #     combine at midpoint
    V0 = (AKs0u-AKs0b).squeeze()/(z0u-z0b)
    
    # second time step
    ZH1 = get_V(vn_list_zh, ds1, plon, plat, pcs, R, surface)
    ZH1[np.isnan(ZH1)] = 0
    dpcs1 = 1/(ZH1[:,0] + ZH1[:,1]) # change in pcs for 1m difference
    
    #     upper variables
    pcs1u = pcs-dpcs1
    pcs1u[pcs0u > S['Cs_r'][-1]] = S['Cs_r'][-1]
    AKs1u = get_V(['AKs',], ds1, plon, plat, pcs1u, R, surface) # diffusivity
    z1u = (pcs1u)*(ZH1[:,0]+ZH1[:,1]) # depth = pcs * full-depth
    
    #     lower variables
    pcs1b = pcs+dpcs0
    pcs1b[pcs1b < S['Cs_r'][0]] = S['Cs_r'][0]
    AKs1b = get_V(['AKs',], ds1, plon, plat, pcs0b, R, surface) # diffusivity
    z1b = (pcs1b)*(ZH1[:,0]+ZH1[:,1]) # depth = pcs * full-depth
    
    #     combine at midpoint
    V1 = (AKs1u-AKs1b).squeeze()/(z1u-z1b)
    
    # average of timesteps
    V = (1-frac)*V0 + frac*V1
    
    return V
    
def get_turb(ds0, ds1, dAKs, delta_t, plon, plat, pcs, R, frac, surface):
    # get the vertical turbulence correction components
    
    # getting diffusivity
    V0 = get_V(['AKs',], ds0, plon, plat, pcs, R, surface).squeeze()
    V1 = get_V(['AKs',], ds1, plon, plat, pcs, R, surface).squeeze()
    # replace nans
    V0[np.isnan(V0)] = 0.0
    V1[np.isnan(V1)] = 0.0
    # create weighted average diffusivity
    Vave = (1 - frac)*V0 + frac*V1
    
    # turbulence calculation from Banas, MacCready, and Hickey (2009)
    # w_turbulence = rand*sqrt(2k/delta(t)) + delta(k)/delta(z)
    # rand = random array with normal distribution
    rand = np.random.standard_normal(len(V0))
    # only using half of gradient, first half already added
    V = rand*np.sqrt(2*Vave/delta_t) + dAKs/2
    
    return V

def change_position(ds, plon, plat, pcs, R, surface):
    # if the current position is on the boundary, find the next deepest
    # position so velocity is not 0
    # function will loop through coordinates in a progressively larger
    # radius until a non-boundary position is found

    v_list = ['u', 'v']
    z_list = ['zeta', 'h']
    plon_new = plon.copy()
    plat_new = plat.copy()
    # limit plon/plat to unique coordinates to improve efficiency
    puni, uni_ind, uni_recon = np.unique([' '.join(map(str,j)) for j in 
                            np.array((plon,plat)).T], return_index=True,
                            return_inverse=True)
    plon_uni = plon[uni_ind]
    plat_uni = plat[uni_ind]
    NP = len(plon_uni)
    
    # get initial velocities and positions of rho grid
    V = get_V(v_list, ds, plon_uni, plat_uni, np.ones(NP), R, surface)
    lat_arr = R['rlatr']
    lon_arr = R['rlonr']
    
    for ind in np.arange(NP):
        # where grid is on boundary, all velocities return as 0
        if (V[ind,:] == 0).any():
        
            # closest index in grid
            lat_ind = np.argsort(abs(lat_arr-plat_uni[ind]))
            lon_ind = np.argsort(abs(lon_arr-plon_uni[ind]))
            # lat/lon in order of closeness
            lat_sort = lat_arr[lat_ind]
            lon_sort = lon_arr[lon_ind]
            
            # run through the neighboring tiles, increasing the square
            # size each time until a non-boundary cell is found
            # i.e. 9, 25, 49... grid points around the original point
            for sqr_size in np.arange(2,int(len(lat_arr)/2)+1):
                # pull this square's grid points
                lat_seg = lat_sort[0:2*sqr_size-1]
                lon_seg = lon_sort[0:2*sqr_size-1]
                lat_grid, lon_grid = np.meshgrid(lat_seg, lon_seg)
                # converted to plon/plat format
                plat_seg = lat_grid.flatten()
                plon_seg = lon_grid.flatten()
                pcs_seg = np.ones((len(plon_seg)))
                
                # find depth and velocities at each location
                ZH = get_V(z_list, ds, plon_seg, plat_seg, pcs_seg, R, surface)
                V_test = get_V(v_list, ds, plon_seg, plat_seg, pcs_seg, R, surface)
                dep_seg = -(ZH[:,0] + ZH[:,1])
                # sort by deepest points
                depth_sort = np.argsort(dep_seg)
                # loop through each point until one is not on the boundary
                # the deepest point is used first because that is the most
                # likely to be non_boundary
                for seg_ind in range(len(plon_seg)):
                    deep_ind = depth_sort[seg_ind]
                    # test new position             
                    if (V_test[deep_ind,:] != 0).all():
                        break
                # set new lat/lon if non-boundary has been found
                if (V_test[deep_ind,:] != 0).all():
                    plat_uni[ind] = plat_seg[deep_ind]
                    plon_uni[ind] = plon_seg[deep_ind]
                    break
                

    # put unique values back into full array
    plat_new = plat_uni[uni_recon]
    plon_new = plon_uni[uni_recon]

    return plon_new, plat_new

# Bartos - end edit

def get_properties(vn_list_other, ds, it, P, plon, plat, pcs, R, surface):
    # find properties at a position
    OTH = get_V(vn_list_other, ds, plon, plat, pcs, R, surface)
    for vn in vn_list_other:
        P[vn][it,:] = OTH[:,vn_list_other.index(vn)]
    this_zeta = OTH[:, vn_list_other.index('zeta')]
    this_h = OTH[:, vn_list_other.index('h')]
    full_depth = this_zeta + this_h
    P['z'][it,:] = pcs * full_depth
    P['zbot'][it,:] = -this_h

    return P

def get_V(vn_list, ds, plon, plat, pcs, R, surface, bound='stop'):

    from warnings import filterwarnings
    filterwarnings('ignore') # skip some warning messages

    # get interpolant arrays
    i0lon_d = dict()
    i1lon_d = dict()
    frlon_d = dict()
    i0lat_d = dict()
    i1lat_d = dict()
    frlat_d = dict()
    for gg in ['r', 'u', 'v']:
        exn = False
        i0lon_d[gg], i1lon_d[gg], frlon_d[gg] = zfun.get_interpolant(
                plon, R['rlon'+gg], extrap_nan=exn)
        i0lat_d[gg], i1lat_d[gg], frlat_d[gg] = zfun.get_interpolant(
                plat, R['rlat'+gg], extrap_nan=exn)
    i0csr, i1csr, frcsr = zfun.get_interpolant(pcs, R['rcsr'], extrap_nan=exn)
    i0csw, i1csw, frcsw = zfun.get_interpolant(pcs, R['rcsw'], extrap_nan=exn)
    NV = len(vn_list)
    NP = len(plon)
    # get interpolated values
    V = np.nan * np.ones((NP,NV))
    vcount = 0
    for vn in vn_list:
        if vn in ['w']:
            i0cs = i0csw
            i1cs = i1csw
            frcs = frcsw
        else:
            i0cs = i0csr
            i1cs = i1csr
            frcs = frcsr
# Bartos - adding 'AKs_turb' and 'dAKs_dz' to list
        if vn in ['salt','temp','zeta','h','Uwind','Vwind','w','AKs','dAKs_dz']:
            gg = 'r'
        elif vn in ['u']:
            gg = 'u'
        elif vn in ['v']:
            gg = 'v'
        i0lat = i0lat_d[gg]
        i1lat = i1lat_d[gg]
        frlat = frlat_d[gg]
        i0lon = i0lon_d[gg]
        i1lon = i1lon_d[gg]
        frlon = frlon_d[gg]
        # get the data field and put nan's in masked points
        if vn in ['salt','temp','u','v','w'] and surface==True:
            v0 = ds[vn][0, -1, :, :].squeeze()
        else:
            v0 = ds[vn][:].squeeze()
        try:
            vv = v0.data
            vv[v0.mask] = np.nan
        except AttributeError:
            # it is not a masked array
            vv = v0
# Bartos - adding 'AKs' to list
        if vn in ['salt','temp','AKs','u','v','w'] and surface==False:
            # For variables with depth axis
            # Get just the values around our particle positions.
            # each row in VV corresponds to a "box" around a point
            VV = np.nan* np.ones((NP, 8))
            VV[:,0] = vv[i0cs, i0lat, i0lon]
            VV[:,1] = vv[i0cs, i0lat, i1lon]
            VV[:,2] = vv[i0cs, i1lat, i0lon]
            VV[:,3] = vv[i0cs, i1lat, i1lon]
            VV[:,4] = vv[i1cs, i0lat, i0lon]
            VV[:,5] = vv[i1cs, i0lat, i1lon]
            VV[:,6] = vv[i1cs, i1lat, i0lon]
            VV[:,7] = vv[i1cs, i1lat, i1lon]
            # Work on edge values.  If all in a box are masked
            # then that row will be nan's, and also:
# Bartos - adding reflective boundary condition and 'AKs_turb' to velocity list
            if bound == 'stop':
                if vn in ['u', 'v', 'w', 'AKs']:
                    # set all velocities to zero if any in the box are masked
                    VV[np.isnan(VV).any(axis=1), :] = 0
            elif bound == 'reflect':
                # get magnitude of average velocity
                newval_mag = abs(np.nanmean(VV, axis=1).reshape(NP,1))
                if vn in ['w', 'AKs']:
                    # set vertical velocity and AKs to zero if any in the box are masked
                    VV[np.isnan(VV).any(axis=1), :] = 0
                elif vn in ['u']:
                    # mask0 for masked points on i0lon
                    mask0 = np.isnan(VV[:,[0,2,4,6]]).any(axis=1)
                    # mask1 for masked points on i1lon
                    mask1 = np.isnan(VV[:,[1,3,5,7]]).any(axis=1)
                    # set u-velocity away from mask
                    VV[mask0, :] = newval_mag[mask0]
                    VV[mask1, :] = -newval_mag[mask1]
                elif vn in ['v']:
                    # mask0 for masked points on i0lat
                    mask0 = np.isnan(VV[:,[0,1,4,5]]).any(axis=1)
                    # mask1 for masked points on i1lat
                    mask1 = np.isnan(VV[:,[2,3,6,7]]).any(axis=1)
                    # set v-velocity away from mask
                    VV[mask0, :] = newval_mag[mask0]
                    VV[mask1, :] = -newval_mag[mask1]
            elif vn in ['salt', 'temp']:
                # set all tracers to their average if any in the box are masked
                newval = np.nanmean(VV, axis=1).reshape(NP, 1) * np.ones((1,8))
                mask = np.isnan(VV)
                VV[mask] = newval[mask]
            # now do the interpolation in each box
            vl = ( (1-frlat)*((1-frlon)*VV[:,0] + frlon*VV[:,1])
                + frlat*((1-frlon)*VV[:,2] + frlon*VV[:,3]) )
            vu = ( (1-frlat)*((1-frlon)*VV[:,4] + frlon*VV[:,5])
                + frlat*((1-frlon)*VV[:,6] + frlon*VV[:,7]) )
            v = (1-frcs)*vl + frcs*vu
        elif vn in ['salt','temp','u','v','w'] and surface==True:
            VV = np.nan* np.ones((NP, 4))
            VV[:,0] = vv[i0lat, i0lon]
            VV[:,1] = vv[i0lat, i1lon]
            VV[:,2] = vv[i1lat, i0lon]
            VV[:,3] = vv[i1lat, i1lon]
            # Work on edge values.  If all in a box are masked
            # then that row will be nan's, and also:
# Bartos - adding reflective boundary condition and 'AKs_turb' to velocity list
            if bound == 'stop':
                if vn in ['u', 'v', 'w', 'AKs']:
                    # set all velocities to zero if any in the box are masked
                    VV[np.isnan(VV).any(axis=1), :] = 0
            elif bound == 'reflect':
                # get magnitude of average velocity
                newval_mag = abs(np.nanmean(VV, axis=1).reshape(NP,1))
                if vn in ['w', 'AKs']:
                    # set vertical velocity and AKs to zero if any in the box are masked
                    VV[np.isnan(VV).any(axis=1), :] = 0
                elif vn in ['u']:
                    # mask0 for masked points on i0lon
                    mask0 = np.isnan(VV[:,[0,2]]).any(axis=1)
                    # mask1 for masked points on i1lon
                    mask1 = np.isnan(VV[:,[1,3]]).any(axis=1)
                    # set u-velocity away from mask
                    VV[mask0, :] = newval_mag[mask0]
                    VV[mask1, :] = -newval_mag[mask1]
                elif vn in ['v']:
                    # mask0 for masked points on i0lat
                    mask0 = np.isnan(VV[:,[0,1]]).any(axis=1)
                    # mask1 for masked points on i1lat
                    mask1 = np.isnan(VV[:,[2,3]]).any(axis=1)
                    # set v-velocity away from mask
                    VV[mask0, :] = newval_mag[mask0]
                    VV[mask1, :] = -newval_mag[mask1]
            elif vn in ['salt','temp']:
                # set all tracers to their average if any in the box are masked
                newval = np.nanmean(VV, axis=1).reshape(NP, 1) * np.ones((1,4))
                mask = np.isnan(VV)
                VV[mask] = newval[mask]
            newval = np.nanmean(VV, axis=1).reshape(NP, 1) * np.ones((1,4))
            mask = np.isnan(VV)
            VV[mask] = newval[mask]
            v = ( (1-frlat)*((1-frlon)*VV[:,0] + frlon*VV[:,1])
                + frlat*((1-frlon)*VV[:,2] + frlon*VV[:,3]) )
        elif vn in ['zeta','Uwind','Vwind', 'h']:
            # For variables without depth axis
            VV = np.nan* np.ones((NP, 4))
            VV[:,0] = vv[i0lat, i0lon]
            VV[:,1] = vv[i0lat, i1lon]
            VV[:,2] = vv[i1lat, i0lon]
            VV[:,3] = vv[i1lat, i1lon]
            newval = np.nanmean(VV, axis=1).reshape(NP, 1) * np.ones((1,4))
            mask = np.isnan(VV)
            VV[mask] = newval[mask]
            v = ( (1-frlat)*((1-frlon)*VV[:,0] + frlon*VV[:,1])
                + frlat*((1-frlon)*VV[:,2] + frlon*VV[:,3]) )
        V[:,vcount] = v
        vcount += 1

    return V

# Bartos - added keywords yr to provide a year for models with multiple years
#           and days_to_track with a default of daily runs
def get_fn_list(idt, Ldir, days_to_track=1, yr=None):
    # make the list of input history files

    # LiveOcean version
    if Ldir['gtagex'] in ['cascadia1_base_lo1', 'cascadia1_base_lobio1']:
        date_list = []
        for nday in range(days_to_track):
            fdt = idt + timedelta(nday)
            date_list.append(fdt.strftime('%Y.%m.%d'))
        fn_list = []
        for dd in date_list:
            indir = (Ldir['roms'] + 'output/' + Ldir['gtagex'] +
                    '/f' + dd + '/')
            for hh in range(2,26):
                hhhh = ('0000' + str(hh))[-4:]
                fn_list.append(indir + 'ocean_his_' + hhhh + '.nc')

    # Other ROMS runs version
    elif Ldir['gtagex'] == 'D2005_his':
        # Must be run on Parker's Mac
        indir = '/Users/PM5/Documents/roms/output/' + Ldir['gtagex'] + '/'
        save_num_list = range(1,365*24)
        save_dt_list = []
        dt00 = datetime(2005,1,1)
        save_dt_list.append(dt00)
        for sn in save_num_list:
            save_dt_list.append(dt00 + timedelta(hours=sn))
        # keys of this dict are datetimes, and values are history numbers
        save_dt_num_dict = dict(zip(save_dt_list,save_num_list))
        fn_list = []
# Bartos - changed days_to_track to kwarg, not Ldir, and removed one hour
        for hh in range(days_to_track*24):
            hh = save_dt_num_dict[idt + timedelta(hours=hh)]
            hhhh = ('0000' + str(hh))[-4:]
            fn_list.append(indir + 'ocean_his_' + hhhh + '.nc')
            
# Bartos - added both C2009 and models for rockfish

    elif Ldir['gtagex'] == 'C2009':
        # Must be run on Fjord
        indir = '/boildat1/parker/roms/output/C2009/OUT/'
        save_num_list = range(1,365*24)
        save_dt_list = []
        dt00 = datetime(2009,1,1)
        save_dt_list.append(dt00)
        for sn in save_num_list:
            save_dt_list.append(dt00 + timedelta(hours=sn))
        # keys of this dict are datetimes, and values are history numbers
        save_dt_num_dict = dict(zip(save_dt_list,save_num_list))
        fn_list = []
        for hh in range(days_to_track*24):
            hh = save_dt_num_dict[idt + timedelta(hours=hh)]
            hhhh = ('0000' + str(hh))[-4:]
            fn_list.append(indir + 'ocean_his_' + hhhh + '.nc')
            
    elif Ldir['ic_name'] == 'rockfish':
        # Must be run on Fjord or Boiler
        which_home = os.path.expanduser('~')
        # experiment grid determines indir
        if Ldir['gtagex'] == 'MoSSea':
            indir = '/pmr3/pmraid1/daves/runs/salish_2006_4/OUT/'
        if Ldir['gtagex'] == 'PNWTOX':
            # determine the current computer
            if 'fjord.txt' in os.listdir(which_home):
                ddir = '/boildat1/'
            elif 'boiler.txt' in os.listdir(which_home):
                ddir = '/data1/'
            else:
                raise FileNotFoundError('Need identifier, either "fjord.txt" or "boiler.txt", in home directory.')
            # 2002 pulled from different directory name
            if yr == 2002:
                indir = ddir + 'parker/roms/output/B2002/OUT/'
            else:
                indir = ddir + 'parker/roms/output/C' + str(yr) + '/OUT/'
        # form dictionary matching dates and numbers
        save_num_list = range(1,365*24)
        save_dt_list = []
        dt00 = datetime(yr,1,1)
        save_dt_list.append(dt00)
        for sn in save_num_list:
            save_dt_list.append(dt00 + timedelta(hours=sn))
        # keys of this dict are year then datetime, and values are history numbers
        save_dt_num_dict = dict(zip(save_dt_list,save_num_list))
        fn_list=[]
        # pull all hours, starting from idt
        for hh in range(days_to_track*24):
            hh = save_dt_num_dict[idt + timedelta(hours=hh)]
            hhhh = ('0000' + str(hh))[-4:]
            fn_list.append(indir + 'ocean_his_' + hhhh + '.nc')

    return fn_list