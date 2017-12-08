"""
Plot saved results from tef_0.py.
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
out_fn = out_dir + 'tef_0.p'
D = pickle.load(open(out_fn, 'rb'))
for vn in D.keys():
    locals()[vn] = D[vn]

# derived quantities
td = (T_arr - T_arr[0])/86400 # time in days
dt = T_arr[1] - T_arr[0] # time step in seconds
# these are all
dV_dt = np.diff(V_arr)/dt
dSalt_dt = np.diff(Salt_arr)/dt

q = -(q0_arr + q1_arr)
qs = -(qs0_arr + qs1_arr)

# River flow (positive)
Qr = -Qin1

dVdt_lp = zfun.filt_godin(dV_dt)
dSaltdt_lp = zfun.filt_godin(dSalt_dt)

# PLOTTING
#plt.close('all')

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

#figsize = (18, 12)
figsize = (14, 8)

if True:
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
    ax.set_ylim(-.6, .9)
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
    ax.set_ylim(0, 20)
    ax.set_title('(c) TEF Volume Transports $(10^3\ m^{3}s^{-1})$')
    ax.set_xlabel('Time (days)')
    
    # add spring-neap labels
    add_sn(ax, y=2.5)
    
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
    
plt.show()

# RC CLEANUP
plt.rcdefaults()
