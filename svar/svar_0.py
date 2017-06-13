"""
Code to explore the salinity variance in a simulation.
"""

import os
import sys
import netCDF4 as nc
import matplotlib.pyplot as plt
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
V_arr = np.nan * np.ones(nt+1)
Salt_arr = np.nan * np.ones(nt+1)
SV_arr = np.nan * np.ones(nt+1)
Mix_arr = np.nan * np.ones(nt+1)
# initialize result arrays for average and diagnostic files
T_arr = np.nan + np.ones(nt)
q0_arr = np.nan * np.ones(nt)
q1_arr = np.nan * np.ones(nt)
qs0_arr = np.nan * np.ones(nt)
qs1_arr = np.nan * np.ones(nt)
qsv0_arr = np.nan * np.ones(nt)
qsv1_arr = np.nan * np.ones(nt)
sbar_arr = np.nan * np.ones(nt)

# initialize intermediate results arrays for TEF quantities
sedges = np.linspace(0, 35, 35*20 + 1)
sbins = sedges[:-1] + np.diff(sedges)/2
ns = len(sbins) # number of salinity bins
tef_q0 = np.zeros((ns, nt))
tef_q1 = np.zeros((ns, nt))
tef_qs0 = np.zeros((ns, nt))
tef_qs1 = np.zeros((ns, nt))
tef_qsv0 = np.zeros((ns, nt))
tef_qsv1 = np.zeros((ns, nt))

# loop over hours
tt = 0
tta = 0
for f_dir in f_list:
    print(str(tt))
    
    h_list = os.listdir(dir0 + f_dir)
    h_list.sort()
    h_list = [x for x in h_list if x[:9]=='ocean_his']
    if tt == 0:
        pass
    else:
        h_list = h_list[1:] # drop the zero hour for all but the first day
        
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
        sbar = Salt/V # average salinity
        sp = salt - sbar # s' = s - sbar
        sv = sp**2 # salinity variance
        SV = np.sum(da.reshape((1,ny,nx))*dzr*sv)
                
        # and calculate net destruction of variance by vertical mixing
        K = ds['AKs'][0, 1:-1, j0:j1, i0:i1].squeeze()
        dzw = np.diff(zr, axis=0)
        dsdz = np.diff(salt, axis=0) / dzw
        mix = -2 * K * dsdz**2
        dvw = da.reshape((1,ny,nx))*dzw
        Mix = np.sum(mix * dvw)
        
        # store results
        V_arr[tt] = V
        Salt_arr[tt] = Salt
        SV_arr[tt] = SV
        Mix_arr[tt] = Mix
        
        ds.close()
        tt += 1

    a_list = os.listdir(dir0 + f_dir)
    a_list.sort()
    a_list = [x for x in a_list if x[:9]=='ocean_avg']
    for ai in a_list:
        
        fn = dir0 + f_dir + '/' + ai
        
        ds = nc.Dataset(fn)
        
        T_arr[tta] = ds['ocean_time'][:].squeeze()
        
        dq0 = ds['Huon'][0, : , j0:j1, i0-1].squeeze()
        dq1 = ds['Huon'][0, : , j0:j1, i1-1].squeeze()
        
        dqs0 = ds['Huon_salt'][0, : , j0:j1, i0-1].squeeze()
        dqs1 = ds['Huon_salt'][0, : , j0:j1, i1-1].squeeze()
        
        # calculating flux of salinity variance
        #
        # first we need to get sbar for the average state
        zeta = ds['zeta'][0, j0:j1, i0:i1].squeeze()
        salt = ds['salt'][0, :, j0:j1, i0:i1].squeeze()
        zr, zw = zrfun.get_z(h, zeta, S)
        dzr = np.diff(zw, axis=0)        
        V = np.sum(da.reshape((1,ny,nx))*dzr) # volume
        Salt = np.sum(da.reshape((1,ny,nx))*dzr*salt) # net salt
        sbar = Salt/V # average salinity
        sbar_arr[tta] = sbar
        # then get the salinity averaged onto the u-grid on both open boundaries
        s0 = (ds['salt'][0, :, j0:j1, i0-1].squeeze() + ds['salt'][0, :, j0:j1, i0].squeeze())/2
        s1 = (ds['salt'][0, :, j0:j1, i1-1].squeeze() + ds['salt'][0, :, j0:j1, i1].squeeze())/2
        # and turn it into variance
        sp0 = s0 - sbar
        sp1 = s1 - sbar
        sv0 = sp0**2
        sv1 = sp1**2

        # Flux calculations
        #
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
        # Salinity Variance
        qsv0_arr[tta] = -np.sum(dq0 * sv0)
        qsv1_arr[tta] = np.sum(dq1 * sv1)
                
        # TEF variables
        s00 = s0[s0.mask==False] # also flattens the array
        dq00 = dq0[dq0.mask==False]
        dqs00 = dqs0[dqs0.mask==False]
        dqsv0 = dq0 * sv0
        dqsv00 = dqsv0[dqsv0.mask==False]
        inds = np.digitize(s00, sedges, right=True)
        counter = 0
        for ii in inds:
            tef_q0[ii-1,tta] += dq00[counter]
            tef_qs0[ii-1,tta] += dqs00[counter]
            tef_qsv0[ii-1,tta] += dqsv00[counter]
            counter += 1
        #
        s11 = s1[s1.mask==False]
        dq11 = dq1[dq1.mask==False]
        dqs11 = dqs1[dqs1.mask==False]
        dqsv1 = dq1 * sv1
        dqsv11 = dqsv1[dqsv1.mask==False]
        inds = np.digitize(s11, sedges, right=True)
        counter = 0
        for ii in inds:
            tef_q1[ii-1,tta] += dq11[counter]
            tef_qs1[ii-1,tta] += dqs11[counter]
            tef_qsv1[ii-1,tta] += dqsv11[counter]
            counter += 1
        ds.close()
        
        tta += 1

#%% TEF processing

tef_q0_lp = np.nan * np.ones_like(tef_q0)
tef_q1_lp = np.nan * np.ones_like(tef_q1)
tef_qs0_lp = np.nan * np.ones_like(tef_qs0)
tef_qs1_lp = np.nan * np.ones_like(tef_qs1)
tef_qsv0_lp = np.nan * np.ones_like(tef_qsv0)
tef_qsv1_lp = np.nan * np.ones_like(tef_qsv1)
for ii in range(ns):
    tef_q0_lp[ii,:] = zfun.filt_godin(tef_q0[ii,:])
    tef_q1_lp[ii,:] = zfun.filt_godin(tef_q1[ii,:])
    tef_qs0_lp[ii,:] = zfun.filt_godin(tef_qs0[ii,:])
    tef_qs1_lp[ii,:] = zfun.filt_godin(tef_qs1[ii,:])
    tef_qsv0_lp[ii,:] = zfun.filt_godin(tef_qsv0[ii,:])
    tef_qsv1_lp[ii,:] = zfun.filt_godin(tef_qsv1[ii,:])
   
qin = tef_q0_lp.copy()
qout = tef_q0_lp.copy()
qsin = tef_qs0_lp.copy()
qsout = tef_qs0_lp.copy()
qsvin = tef_qsv0_lp.copy()
qsvout = tef_qsv0_lp.copy()

qin[tef_q0_lp<0] = 0
qout[tef_q0_lp>0] = 0
qsin[tef_q0_lp<0] = 0
qsout[tef_q0_lp>0] = 0
qsvin[tef_q0_lp<0] = 0
qsvout[tef_q0_lp>0] = 0
#
Qin0 = qin.sum(axis=0)
Qout0 = qout.sum(axis=0)
QSin0 = qsin.sum(axis=0)
QSout0 = qsout.sum(axis=0)
QSVin0 = qsvin.sum(axis=0)
QSVout0 = qsvout.sum(axis=0)
#
Sin0 = QSin0/Qin0
Sout0 = QSout0/Qout0

SVin0 = QSVin0/Qin0
SVout0 = QSVout0/Qout0

qin = tef_q1_lp.copy()
qout = tef_q1_lp.copy()
qsin = tef_qs1_lp.copy()
qsout = tef_qs1_lp.copy()
qsvin = tef_qsv1_lp.copy()
qsvout = tef_qsv1_lp.copy()
#
qin[tef_q1_lp>0] = 0 # switch signs compared to open boundary 0
qout[tef_q1_lp<0] = 0
qsin[tef_q1_lp>0] = 0
qsout[tef_q1_lp<0] = 0
qsvin[tef_q1_lp>0] = 0
qsvout[tef_q1_lp<0] = 0
#
Qin1 = qin.sum(axis=0)
Qout1 = qout.sum(axis=0)
QSin1 = qsin.sum(axis=0)
QSout1 = qsout.sum(axis=0)
QSVin1 = qsvin.sum(axis=0)
QSVout1 = qsvout.sum(axis=0)
#
#
Sin1 = QSin1/Qin1
#Sout1 = QSout1/Qout1 # Qout1 = 0
SVin1 = QSVin1/Qin1
#SVout1 = QSVout1/Qout1

#%% Plotting

td = (T_arr - T_arr[0])/86400 # time in days
dt = T_arr[1] - T_arr[0] # time step in seconds

dV_dt = np.diff(V_arr)/dt
dSalt_dt = np.diff(Salt_arr)/dt
dSV_dt = np.diff(SV_arr)/dt

q = -(q0_arr + q1_arr)
qs = -(qs0_arr + qs1_arr)
qsv = -(qsv0_arr + qsv1_arr)

Mixa_arr = Mix_arr[:-1] + np.diff(Mix_arr)/2

Sbar = zfun.filt_godin(sbar_arr)

V0 = V_arr.mean() # mean volume
SV_max = (35/2)**2 * V0 # max volume-integrated salinity variance
dSV_dt_scale = SV_max / (7*86400)

# River flow (positive)
Qr = -Qin1

td0 = td[0]
td1 = td[-1]

plt.close('all')

if True:
    # TEF
    fig0 = plt.figure(figsize=(13,8))
    
    ax = fig0.add_subplot(2,2,1)
    Sstorage = zfun.filt_godin(dSalt_dt)
    ax.plot(td, Sstorage/dSV_dt_scale, '-r', label='dS/dt')
    ax.plot(td, Qin0*Sin0/dSV_dt_scale, '-b', label='Qin*Sin')
    ax.plot(td, Qout0*Sout0/dSV_dt_scale, '--r', label='Qout*Sout')
    ax.plot(td, (Sstorage - Qin0*Sin0 - Qout0*Sout0)/dSV_dt_scale, '-k', label='Error')
    ax.legend()
    ax.grid()
    ax.set_xlim(td0, td1)
    ax.set_title('(a) TEF Salinity Budget')
    
    ax = fig0.add_subplot(2,2,2)
    ax.plot(td, q0_arr/1000, '-m')
    ax.set_xlim(td0, td1)
    ax.set_title('(b) Tidal Volume Transport at Mouth $(1000\ m^{3}s^{-1})$')
    
    ax = fig0.add_subplot(2,2,3)
    ax.plot(td, Qin0, '-r')
    ax.plot(td, -Qout0, '-b')
    ax.plot(td, -Qin1, '--r')
    #ax.plot(td, Qout1, '--b')
    ax.legend(['Qin','-Qout','Qr'])
    ax.grid()
    ax.set_xlim(td0, td1)
    ax.set_title('(c) TEF Volume Transports')
    ax.set_xlabel('Time (days)')
    
    ax = fig0.add_subplot(2,2,4)
    ax.plot(td, Sin0, '-r')
    ax.plot(td, Sout0, '-b')
    ax.legend(['Sin','Sout'])
    ax.grid()
    ax.set_xlim(td0, td1)
    ax.set_title('(d) TEF Salinities')
    ax.set_xlabel('Time (days)')
    
if False:
    # plot raw time series of volume and salt equations
    fig1 = plt.figure(figsize=(13,8))
    #
    ax = fig1.add_subplot(2,1,1)
    ax.plot(td, q, '-b')
    ax.plot(td, dV_dt, '-r')
    ax.set_title('Test of Volume Conservation')
    #
    ax = fig1.add_subplot(2,1,2)
    ax.plot(td, qs, '-b')
    ax.plot(td, dSalt_dt, '-r')
    ax.set_title('Test of Salt Conservation')
    ax.set_xlabel('Time (days)')
    #
    rmse = np.std(qs-dSalt_dt)
    relative_error = rmse / np.std(dSalt_dt)
    print('rms salt error / std(dSalt/dt) = %0.5f' % (relative_error))
    
if True:
    # plot low-passed time series of salinity variance equation
    fig2 = plt.figure(figsize=(13,8))
    
    ax = fig2.add_subplot(2,1,1)
    l1, = ax.plot(td, zfun.filt_godin(dSV_dt)/dSV_dt_scale, '-r')
    l1.set_label('dSV/dt')
    l2, = ax.plot(td, zfun.filt_godin(qsv)/dSV_dt_scale, '-b')
    l2.set_label('Advection')
    l3, = ax.plot(td, zfun.filt_godin(Mixa_arr)/dSV_dt_scale, 'g')
    l3.set_label('Mixing')
    l4, = ax.plot(td, zfun.filt_godin(dSV_dt - qsv - Mixa_arr)/dSV_dt_scale, '-k')
    l4.set_label('Error = Numerical Mixing')
    ax.legend()
    ax.set_title('(a) Low-passed Salt Variance: dSV/dt = Advection + Mixing + Error')
    ax.grid()
    ax.set_ylabel('Terms $(SV_{max}\ (7 days)^{-1}))$')
    yl = .65
    ax.set_ylim(-yl, yl)
    ax.set_xlim(td0, td1)    
    
    ax = fig2.add_subplot(2,1,2)    
    l5, = ax.plot(td, Qin0*SVin0/dSV_dt_scale, '-m')
    l5.set_label('Qin*SVin')
    l6, = ax.plot(td, Qout0*SVout0/dSV_dt_scale, ':m')
    l6.set_label('Qout*SVout')
    l7, = ax.plot(td, Qr*SVin1/dSV_dt_scale, '--m')
    l7.set_label('Qr*SVin at river end')
    l8, = ax.plot(td, (Qin0*SVin0+Qout0*SVout0+Qr*SVin1)/dSV_dt_scale, '-b')
    l8.set_label('Sum = Advection')
    
    ax.legend()
    ax.set_title('(b) Advection Decomposed into TEF Terms')
    ax.grid()
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Terms $(SV_{max}\ (7 days)^{-1}))$')
    ax.set_ylim(-yl, yl)
    ax.set_xlim(td0, td1)
    
if True:
    # plot low-passed time series of salinity variance equation
    fig3 = plt.figure(figsize=(13,8))
    
    yl = .65

    ax = fig3.add_subplot(1,1,1)    

    l8, = ax.plot(td, (Qin0*SVin0+Qout0*SVout0+Qr*SVin1)/dSV_dt_scale, '-b')
    l8.set_label('Correct Advection')
    
    alt1 = Sin0*Sout0*Qr
    alt2 = (Sin0+Sout0-2*Sbar)*Sstorage
    
    l9, = ax.plot(td, alt1/dSV_dt_scale, ':g')
    l9.set_label('Sin*Sout*Qr')
    l10, = ax.plot(td, alt2/dSV_dt_scale, '--g')
    l10.set_label('(Sin+Sout-2Sbar)*dS/dt')
    l11, = ax.plot(td, (alt1 + alt2)/dSV_dt_scale, '-g')
    l11.set_label('Sum = Approximate Advection')
    
    ax.legend()
    ax.set_title('Approximate Advection Terms')
    ax.grid()
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Terms $(SV_{max}\ (7 days)^{-1}))$')
    ax.set_ylim(-yl, yl)
    ax.set_xlim(td0, td1)
    
plt.show()

