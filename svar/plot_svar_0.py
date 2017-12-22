"""
Plot saved results from svar_0.py.
"""

import pickle
import matplotlib.pyplot as plt
import numpy as np

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import zfun
import Lfun
Ldir = Lfun.Lstart()

# get saved fields
out_dir = Ldir['parent'] + 'ptools_output/svar/'
Lfun.make_dir(out_dir)
out_fn = out_dir + 'svar_0.p'
D = pickle.load(open(out_fn, 'rb'))
for vn in D.keys():
    locals()[vn] = D[vn]

# derived quantities
td = (T_arr - T_arr[0])/86400 # time in days
dt = T_arr[1] - T_arr[0] # time step in seconds
# these are all
dV_dt = np.diff(V_arr)/dt
dSalt_dt = np.diff(Salt_arr)/dt
dSV_dt = np.diff(SV_arr)/dt

q = -(q0_arr + q1_arr)
qs = -(qs0_arr + qs1_arr)
qsv = -(qsv0_arr + qsv1_arr)

Mixa_arr = Mix_arr[:-1] + np.diff(Mix_arr)/2

Sbar = zfun.filt_godin(sbar_arr)
SVh = zfun.filt_godin(SV_arr)
SV = SVh[:-1] + np.diff(SVh)/2
Vh = zfun.filt_godin(V_arr)
V = Vh[:-1] + np.diff(Vh)/2

V0 = V_arr.mean() # mean volume

# River flow (positive)
Qr = -Qin1

dVdt_lp = zfun.filt_godin(dV_dt)
dSaltdt_lp = zfun.filt_godin(dSalt_dt)
dSVdt_lp = zfun.filt_godin(dSV_dt)
Mix_lp = zfun.filt_godin(Mixa_arr)
MixFull_lp = -zfun.filt_godin(dSV_dt - qsv)
MixTEF_lp = (Sin0*Sout0*Qr +
        dSaltdt_lp * (Sin0 + Sout0 - 2*Sbar) +
        dVdt_lp * (Sbar*Sbar - Sin0*Sout0) -
        dSVdt_lp)
MixTEFsimple = Sin0*Sout0*Qr

AdvFull = dSVdt_lp + MixFull_lp
AdvTEF = dSVdt_lp + MixTEF_lp


# record-mean budgets
MM = np.nanmean(zfun.filt_godin(Mix_arr))
MM_full = np.nanmean(MixFull_lp)
SSin = np.nanmean(Qin0*Sin0)/np.nanmean(Qin0)
SSout = np.nanmean(Qout0*Sout0)/np.nanmean(Qout0)
print('SSin = %0.2f and SSout = %0.2f' % (SSin, SSout))
QQr = np.nanmean(Qr)
MM_alt = QQr * SSin * SSout
print('MM_full = %0.2f, Qr*Sin*Sout = %0.2f (1e6 (g/kg)^2 m3/s)' % (MM_full/1e6, MM_alt/1e6))

# PLOTTING
plt.close('all')

# time limits for plotting
td0 = td[0]
td1 = td[-1]

# RC SETUP (plotting defaults)
def set_rc(fs_big, fs_small, lw_big, lw_small):
    plt.rc('xtick', labelsize=fs_small)
    plt.rc('ytick', labelsize=fs_small)
    plt.rc('xtick.major', size=10, pad=5, width=lw_small)
    plt.rc('ytick.major', size=10, pad=5, width=lw_small)
    plt.rc('axes', lw=lw_small)
    plt.rc('lines', lw=lw_big)
    plt.rc('font', size=fs_big)
    plt.rc('grid', color='g', ls='-', lw=lw_small, alpha=.3)
fs_big = 12
fs_small = 10
lw_big = 3
lw_small = 2
set_rc(fs_big, fs_small, lw_big, lw_small)

# a function to add spring neap labels
def add_sn(ax, y=0):
    # add spring-neap labels
    snx_list = [3, 10.4, 17.8, 25.2]
    sn_list = ['S', 'N', 'S', 'N']
    sn_dict = dict(zip(snx_list, sn_list))
    for snx in sn_dict.keys():
        sn = sn_dict[snx]
        ax.text(snx, y, sn, fontsize=50, alpha=.5, horizontalalignment='center')

figsize = (14, 12)

if False:
    # TEF Salt Budget
    fig0 = plt.figure(figsize=figsize)
    
    ax = fig0.add_subplot(2,2,1)
    scl = 1e6
    Sstorage = zfun.filt_godin(dSalt_dt)
    ax.plot(td, Sstorage/scl, '-r', label='d(Net Salt)/dt')
    ax.plot(td, Qin0*Sin0/scl, '-m', label='$Q_{in}S_{in}$')
    ax.plot(td, Qout0*Sout0/scl, '--m', label='$Q_{out}S_{out}$')
    ax.plot(td, (Sstorage - Qin0*Sin0 - Qout0*Sout0)/scl, '-k', label='Error')
    ax.legend(ncol=2, loc='upper center')
    ax.grid()
    ax.set_xlim(td0, td1)
    ax.set_ylim(-.4, .6)
    ax.set_title('(a) TEF Salinity Budget $(10^6\ (g/kg)\ m^{3}s^{-1})$')
    ax.set_xticklabels('')
    ax.set_xlabel('')
    
    ax = fig0.add_subplot(2,2,2)
    scl = 1e3
    ax.plot(td, q0_arr/scl, '-k', linewidth=1.5)
    ax.set_xlim(td0, td1)
    ax.set_ylim(-250, 250)
    ax.grid()
    ax.set_title('(b) Tidal Volume Transport at Mouth $(10^3\ m^{3}s^{-1})$')
    ax.set_xticklabels('')
    ax.set_xlabel('')
    
    # add spring-neap labels
    add_sn(ax, y=-200)
    
    ax = fig0.add_subplot(2,2,3)
    scl = 1e3
    ax.plot(td, Qin0/scl, '-b')
    ax.plot(td, -Qout0/scl, '--b')
    ax.plot(td, -Qin1/scl, ':b')
    ax.legend(['$Q_{in}$','$-Q_{out}$','$Q_r$'], ncol=3, loc='upper center')
    ax.grid()
    ax.set_xlim(td0, td1)
    ax.set_ylim(0, 16)
    ax.set_title('(c) TEF Volume Transports $(10^3\ m^{3}s^{-1})$')
    ax.set_xlabel('Time (days)')
    
    # add spring-neap labels
    add_sn(ax, y=5)
    
    ax = fig0.add_subplot(2,2,4)
    ax.plot(td, Sin0, '-', color='orange')
    ax.plot(td, Sout0, '--', color='orange')
    ax.legend(['$S_{in}$','$S_{out}$'], loc='lower left')
    ax.grid()
    ax.set_xlim(td0, td1)
    ax.set_ylim(20, 35)
    ax.set_title('(d) TEF Salinities')
    ax.set_xlabel('Time (days)')
    # add spring-neap labels
    add_sn(ax, y=26)
    
    
if True:
    # Variance Budget
    fig2 = plt.figure(figsize=figsize)
    
    scl = 1e6
    
    ax = fig2.add_subplot(3,1,1)
    l1, = ax.plot(td, zfun.filt_godin(dSV_dt)/scl, '-r')
    l1.set_label('d(Net Variance)/dt')
    l2, = ax.plot(td, zfun.filt_godin(qsv)/scl, '-b')
    l2.set_label('Advection')
    l3, = ax.plot(td, -zfun.filt_godin(Mixa_arr)/scl, 'g')
    l3.set_label('-Mixing')
    l4, = ax.plot(td, zfun.filt_godin(dSV_dt - qsv + Mixa_arr)/scl, '-k')
    l4.set_label('-Numerical Mixing')
    ax.legend(ncol=4, loc='upper center')
    ax.set_title('(a) d(Net Variance)/dt = Advection - Mixing - Numerical Mixing')
    ax.grid()
    ylab = '$10^6\ (g/kg)^2\ m^3s^{-1}$'
    ax.set_ylabel(ylab)
    ax.set_ylim(-3, 4)
    ax.set_xlim(td0, td1)    
    ax.set_xticklabels('')
    ax.set_xlabel('')
    ax.plot([td0, td1], [0,0], '-k', linewidth=1)
    
    ax = fig2.add_subplot(3,1,2)    
    l5, = ax.plot(td, Qin0*SVin0/scl, '-m')
    l5.set_label('$Q_{in}S^{\prime2}_{in}$')
    l6, = ax.plot(td, Qout0*SVout0/scl, '--m')
    l6.set_label('$Q_{out}S^{\prime2}_{out}$')
    l7, = ax.plot(td, Qr*SVin1/scl, ':m')
    l7.set_label('$Q_{r}\overline{S}^{2}$')
    l8, = ax.plot(td, (Qin0*SVin0+Qout0*SVout0+Qr*SVin1)/scl, '-b')
    l8.set_label('Sum = Advection')
    
    ax.legend(ncol=4, loc='upper center')
    ax.set_title('(b) Advection Decomposed into TEF Terms')
    ax.grid()
    ax.set_ylabel(ylab)
    ax.set_ylim(-2, 4)
    ax.set_xlim(td0, td1)
    ax.set_xticklabels('')
    ax.set_xlabel('')
    ax.plot([td0, td1], [0,0], '-k', linewidth=1)
    
    # Just look at total variance over time
    ax = fig2.add_subplot(3,1,3)
    ax.plot(td, SV/V, '-k')
    ax.grid()
    ax.set_title('(c) Volume-Averaged Variance')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('$(g/kg)^2$')
    ax.set_xlim(td0, td1)
    ax.set_xlabel('Time (days)')
    ax.set_ylim(0, 150)
    
    # add spring-neap labels
    add_sn(ax, y=20)
    
if False:
    # three estimates of mixing
    fig4 = plt.figure(figsize=figsize)
    scl = 1e6
    
    ax = fig4.add_subplot(3,1,1)
    l3, = ax.plot(td, MixTEF_lp/scl, '-c')
    l3.set_label('Approximate Mixing')
    l2, = ax.plot(td,MixFull_lp/scl, '-m')
    l2.set_label('Full Mixing')
    l1, = ax.plot(td, Mix_lp/scl, '-g')
    l1.set_label('Resolved Mixing')
    # l4, = ax.plot(td, MixTEFsimple/scl, '-', color='orange')
    # l4.set_label('$Q_rS_{in}S_{out}$')
    ax.legend(loc='upper left')
    ax.text(.6, .9, '(a) Estimates of Mixing', transform=ax.transAxes)
    ax.grid()
    ylab = '$10^6\ (g/kg)^2\ m^3s^{-1}$'
    ax.set_ylabel(ylab)
    ax.set_ylim(0, 4)
    ax.set_xlim(td0, td1)    
    ax.set_xticklabels('')
    ax.set_xlabel('')
    ax.plot([td0, td1], [0,0], '-k', linewidth=1)
    
    # add spring-neap labels
    add_sn(ax, y=.5)
    
    ax = fig4.add_subplot(3,1,2)
    l2, = ax.plot(td, AdvTEF/scl, '-c')
    l2.set_label('Approximate Advection')
    l1, = ax.plot(td, AdvFull/scl, '-b')
    l1.set_label('Advection')
    ax.legend(loc='upper left')
    ax.text(.6, .9, '(b) Estimates of Advection', transform=ax.transAxes)
    ax.grid()
    ylab = '$10^6\ (g/kg)^2\ m^3s^{-1}$'
    ax.set_ylabel(ylab)
    ax.set_ylim(0, 4)
    ax.set_xlim(td0, td1)    
    ax.set_xticklabels('')
    ax.set_xlabel('')
    ax.plot([td0, td1], [0,0], '-k', linewidth=1)
    #
    # compare ways of calculating the TEF Salinity Variance terms
    ax = fig4.add_subplot(3,1,3)
    SVin_alt = (Sin0 - Sbar)**2
    SVout_alt = (Sout0 - Sbar)**2
    ax.plot(td, SVin0,'-m', td, SVout0, '--m',
        td, SVin_alt,'-c', td, SVout_alt, '--c',
        linewidth=3)
    ax.grid()
    ax.text(.05, .9, '(c) TEF Variance at Ocean End', transform=ax.transAxes)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('$(g/kg)^2$')
    ax.set_xlim(td0, td1)
    ax.set_ylim(-10, 350)
    ax.legend(['$S^{\prime2}_{in}$', '$S^{\prime2}_{out}$',
        '$(S_{in}-\overline{S})^2$','$(S_{out}-\overline{S})^2$'],
        loc='lower left', ncol=2)
        
plt.show()

# RC CLEANUP
plt.rcdefaults()
