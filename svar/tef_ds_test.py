"""
Code to explore ways to calcualte TEF, motivated by
the lack of convergence Elizabeth found unsing smaller
salinity bins.
"""

import os
import sys
import netCDF4 as nc
import pickle
import numpy as np

pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import zfun
import zrfun
import Lfun
Ldir = Lfun.Lstart('aestus1', 'A1')
Ldir['gtagex'] = Ldir['gtag'] + '_' + 'ae1'

dir0 = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[0]=='f']

testing = False
if testing:
    f_list = f_list[:14] # control number of days

# number of time steps for averages and diagnostics
nt = len(f_list) * 24

# initialize result arrays for history files
# (one more than the number of ave and dia files)
# these are all volume integrals with tidal variability
V_arr = np.nan * np.ones(nt+1)
Salt_arr = np.nan * np.ones(nt+1)
# initialize result arrays for average files
# 0 and 1 mean ocean and river ends
T_arr = np.nan + np.ones(nt) # ocean time
# volume flux
q0_arr = np.nan * np.ones(nt)
q1_arr = np.nan * np.ones(nt)
# salt flux
qs0_arr = np.nan * np.ones(nt)
qs1_arr = np.nan * np.ones(nt)
# volume-average salinity
sbar_arr = np.nan * np.ones(nt)

# initialize intermediate results arrays for TEF quantities
sedges = np.linspace(0, 35, 1001) # original was 35*20 + 1 bins
sbins = sedges[:-1] + np.diff(sedges)/2
ns = len(sbins) # number of salinity bins
tef_q0 = np.zeros((ns, nt))
tef_q1 = np.zeros((ns, nt))
tef_qs0 = np.zeros((ns, nt))
tef_qs1 = np.zeros((ns, nt))

# loop over hours
tt = 0 # history
tta = 0 # averages
for f_dir in f_list:
    print(str(tt))
    
    # get volume-integrated quantities from the history files
    h_list = os.listdir(dir0 + f_dir)
    h_list.sort()
    h_list = [x for x in h_list if x[:9]=='ocean_his']
    if tt == 0:
        pass
    else:
        h_list = h_list[1:] # drop the zero hour for all but the first day
    #
    for hi in h_list:
        fn = dir0 + f_dir + '/' + hi
        ds = nc.Dataset(fn)
        if tt == 0:
            G, S, T = zrfun.get_basic_info(fn)
            lon_vec = G['lon_rho'][0,:].squeeze()
            lat_vec = G['lat_rho'][:,0].squeeze()
            ii0, ii1, ifr = zfun.get_interpolant(np.array([0.02, 1.5]), lon_vec)
            jj0, jj1, jfr = zfun.get_interpolant(np.array([44.9, 45.1]), lat_vec)
            i0 = ii0[0]
            i1 = ii1[1]
            j0 = jj0[0]
            j1 = jj1[1]
            h = G['h'][j0:j1, i0:i1]
            dx = G['DX'][j0:j1, i0:i1]
            dy = G['DY'][j0:j1, i0:i1]
            da = dx*dy
            ny, nx = da.shape
                        
        zeta = ds['zeta'][0, j0:j1, i0:i1].squeeze()
        salt = ds['salt'][0, :, j0:j1, i0:i1].squeeze()
        zr, zw = zrfun.get_z(h, zeta, S)
        dzr = np.diff(zw, axis=0)
        
        V = np.sum(da.reshape((1,ny,nx))*dzr) # volume
        Salt = np.sum(da.reshape((1,ny,nx))*dzr*salt) # net salt
        
        # store results
        V_arr[tt] = V
        Salt_arr[tt] = Salt
        
        ds.close()
        tt += 1
        
    # next get arrays for flux claculations from the averages
    a_list = os.listdir(dir0 + f_dir)
    a_list.sort()
    a_list = [x for x in a_list if x[:9]=='ocean_avg']
    for ai in a_list:
        fn = dir0 + f_dir + '/' + ai
        ds = nc.Dataset(fn)
        T_arr[tta] = ds['ocean_time'][:].squeeze()
        # getting fluxes of volume and salt
        dq0 = ds['Huon'][0, : , j0:j1, i0-1].squeeze()
        dq1 = ds['Huon'][0, : , j0:j1, i1-1].squeeze()
        dqs0 = ds['Huon_salt'][0, : , j0:j1, i0-1].squeeze()
        dqs1 = ds['Huon_salt'][0, : , j0:j1, i1-1].squeeze()
        # then get the salinity averaged onto the u-grid on both open boundaries
        s0 = (ds['salt'][0, :, j0:j1, i0-1].squeeze() + ds['salt'][0, :, j0:j1, i0].squeeze())/2
        s1 = (ds['salt'][0, :, j0:j1, i1-1].squeeze() + ds['salt'][0, :, j0:j1, i1].squeeze())/2
        # next we take area integrals at ocean and river ends to get net fluxes
        # Volume
        q0_arr[tta] = -np.sum(dq0)
        q1_arr[tta] = np.sum(dq1)
        # Salt
        if True:
            # using the averages
            qs0_arr[tta] = -np.sum(dqs0)
            qs1_arr[tta] = np.sum(dqs1)
        else:
            # alternate version, using volume * salinity
            qs0_arr[tta] = -np.sum(dq0 * s0)
            qs1_arr[tta] = np.sum(dq1 * s1)
        # RESULT these two methods gave nearly identical results
        # with the rms error being about 0.2% of the std of dSalt/dt.
        # This gives support to the calculation of variance advection
        # below.

        # TEF variables
        # which are also area integrals at ocean and river ends
        s00 = s0[s0.mask==False] # flattens the array
        dq00 = dq0[dq0.mask==False]
        dqs00 = dqs0[dqs0.mask==False]
        inds = np.digitize(s00, sedges, right=True)
        counter = 0
        for ii in inds:
            tef_q0[ii-1,tta] += dq00[counter]
            tef_qs0[ii-1,tta] += dqs00[counter]
            counter += 1
        #
        s11 = s1[s1.mask==False]
        dq11 = dq1[dq1.mask==False]
        dqs11 = dqs1[dqs1.mask==False]
        inds = np.digitize(s11, sedges, right=True)
        counter = 0
        for ii in inds:
            # at each time step these are vectors of hourly transport in
            # salinity bins (centered at sbins)
            tef_q1[ii-1,tta] += dq11[counter]
            tef_qs1[ii-1,tta] += dqs11[counter]
            counter += 1
        ds.close()
        
        tta += 1

#%% TEF processing
# first form tidal averages
tef_q0_lp = np.nan * np.ones_like(tef_q0)
tef_q1_lp = np.nan * np.ones_like(tef_q1)
tef_qs0_lp = np.nan * np.ones_like(tef_qs0)
tef_qs1_lp = np.nan * np.ones_like(tef_qs1)
for ii in range(ns):
    tef_q0_lp[ii,:] = zfun.filt_godin(tef_q0[ii,:])
    tef_q1_lp[ii,:] = zfun.filt_godin(tef_q1[ii,:])
    tef_qs0_lp[ii,:] = zfun.filt_godin(tef_qs0[ii,:])
    tef_qs1_lp[ii,:] = zfun.filt_godin(tef_qs1[ii,:])
if False:
    qin = tef_q0_lp.copy()
    qout = tef_q0_lp.copy()
    qsin = tef_qs0_lp.copy()
    qsout = tef_qs0_lp.copy()
    # then mask for in and out parts (ocean end)
    qin[tef_q0_lp<0] = 0
    qout[tef_q0_lp>0] = 0
    qsin[tef_q0_lp<0] = 0
    qsout[tef_q0_lp>0] = 0
    # and integrate over salinity
    Qin0 = qin.sum(axis=0)
    Qout0 = qout.sum(axis=0)
    QSin0 = qsin.sum(axis=0)
    QSout0 = qsout.sum(axis=0)
else:
    # alternate method using cumulative sum of the transport
    # to identify the salinity dividing inflow and outflow
    # RESULT: this way is not sensitive to the number of
    # salinity bins.
    rq0 = np.flipud(tef_q0_lp)
    rqs0 = np.flipud(tef_qs0_lp)
    qcs = np.cumsum(rq0, axis=0)
    imax = np.argmax(qcs, axis=0)
    Qin0 = np.zeros(nt)
    QSin0 = np.zeros(nt)
    Qout0 = np.zeros(nt)
    QSout0 = np.zeros(nt)
    for tt in range(nt):
        Qin0[tt] = np.sum(rq0[:imax[tt], tt])
        Qout0[tt] = np.sum(rq0[imax[tt]:, tt])
        QSin0[tt] = np.sum(rqs0[:imax[tt], tt])
        QSout0[tt] = np.sum(rqs0[imax[tt]:, tt])
    # then fix masking
    nmask = np.isnan(tef_q0_lp[0,:])
    Qin0[nmask] = np.nan
    QSin0[nmask] = np.nan
    Qout0[nmask] = np.nan
    QSout0[nmask] = np.nan

# form derived quantities
Sin0 = QSin0/Qin0
Sout0 = QSout0/Qout0

# same steps for the river end
qin = tef_q1_lp.copy()
qout = tef_q1_lp.copy()
qsin = tef_qs1_lp.copy()
qsout = tef_qs1_lp.copy()
#
qin[tef_q1_lp>0] = 0 # switch signs compared to open boundary 0
qout[tef_q1_lp<0] = 0
qsin[tef_q1_lp>0] = 0
qsout[tef_q1_lp<0] = 0
#
Qin1 = qin.sum(axis=0)
Qout1 = qout.sum(axis=0)
QSin1 = qsin.sum(axis=0)
QSout1 = qsout.sum(axis=0)
#
Sin1 = QSin1/Qin1
#Sout1 = QSout1/Qout1 # Qout1 = 0

# Save results for plotting
D = dict()
D_list = ['T_arr', 'V_arr', 'Salt_arr',
    'q0_arr', 'qs0_arr',
    'q1_arr', 'qs1_arr',
    'Qin0', 'Qout0', 'Sin0', 'Sout0',
    'Qin1', 'Sin1']
for vn in D_list:
    D[vn] = locals()[vn]

out_dir = Ldir['parent'] + 'ptools_output/svar/'
Lfun.make_dir(out_dir)
out_fn = out_dir + 'tef_0.p'
pickle.dump(D, open(out_fn, 'wb'))

