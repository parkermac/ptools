"""
Sensitivity analysis for the basin strait model.
"""

import numpy as np
import matplotlib.pyplot as plt
import basin_strait_functions as bsf
from importlib import reload
reload(bsf)

dims, params = bsf.get_dims()
Socn = params['Socn']

landward_basin = True

dt = 86400 # time step (s)
ND = 1000 # number of days
NT = int(ND*dt / 86400)

# forcing
ut = 1
Qr = [100, 200, 500, 1000, 2000, 5000, 10000]
S1q = []; S2q = []; DSq = []; Qinq = []
for qr in Qr:
    # initial condition
    s1 = Socn-2; s2 = Socn-1; s3 = Socn-3; s4 = Socn-2
    # make an equilibrated state
    for nt in range(NT-1):
        s1, s2, s3, s4, qin, qinp = bsf.advance_s(
                qr, ut, dims, params, s1, s2, s3, s4, dt,
                landward_basin=landward_basin)
    ds = s2 - s1
    # store results
    S1q.append(s1); S2q.append(s2); DSq.append(ds); Qinq.append(qin)

# forcing
Ut = [.1, .2, .5, 1, 2, 5, 10]
qr = 1000
S1u = []; S2u = []; DSu = []; Qinu = []
for ut in Ut:
    # initial condition
    s1 = Socn-2; s2 = Socn-1; s3 = Socn-3; s4 = Socn-2
    # make an equilibrated state
    for nt in range(NT-1):
        s1, s2, s3, s4, qin, qinp = bsf.advance_s(
                qr, ut, dims, params, s1, s2, s3, s4, dt,
                landward_basin=landward_basin)
    ds = s2 - s1
    # store results
    S1u.append(s1); S2u.append(s2); DSu.append(ds); Qinu.append(qin)
    
Qr = np.array(Qr)
Ut = np.array(Ut)

S1q = np.array(S1q)
S2q = np.array(S2q)
DSq = np.array(DSq)
Qinq = np.array(Qinq)

S1u = np.array(S1u)
S2u = np.array(S2u)
DSu = np.array(DSu)
Qinu = np.array(Qinu)

# PLOTTING

fs_big = 16
fs_small = 12
lw_big = 3
lw_small = 2
bsf.set_rc(fs_big, fs_small, lw_big, lw_small)

#plt.close('all')
fig = plt.figure(figsize=(12,8))

ax = fig.add_subplot(321)
ax.plot(Qr/1000, S1q, '-b', Qr/1000, S2q, '-r')
ax.axis([0,10,32,20])
ax.text(.05, .8, '$S_{upper}$', horizontalalignment='left',
    transform=ax.transAxes, color='b', fontsize=20)
ax.text(.05, .6, '$S_{lower}$', horizontalalignment='left',
    transform=ax.transAxes, color='r', fontsize=20)
ax.set_xticklabels('')
ax.set_xlabel('')
ax.grid()

ax = fig.add_subplot(323)
ax.plot(Qr/1000, DSq, '-', color='orange')
#
alpha = DSq[-1] / (Qr[-1]**(2/3))
ax.plot(Qr/1000, (alpha*Qr**(2/3)), ':', color='orange')
#
ax.axis([0,10,0,4])
ax.text(.05, .8, '$\Delta S$', horizontalalignment='left',
    transform=ax.transAxes, color='orange', fontsize=20)
ax.set_xticklabels('')
ax.set_xlabel('')
ax.grid()

ax = fig.add_subplot(325)
ax.plot(Qr/1000, Qinq/1000, '-m')
#
alpha = Qinq[-1] / (Qr[-1]**(1/3))
ax.plot(Qr/1000, (alpha*Qr**(1/3))/1000, ':m')
#
ax.axis([0,10,0,70])
ax.text(.05, .8, '$Q_{in}\ (1000\ m^{3}\ s^{-1})$', horizontalalignment='left',
    transform=ax.transAxes, color='m', fontsize=20)
ax.set_xlabel('$Q_{r}\ (1000\ m^{3}\ s^{-1})$')
ax.grid()

ax = fig.add_subplot(322)
ax.plot(Ut, S1u, '-b', Ut, S2u, '-r')
ax.axis([0,10,32,20])
ax.set_xticklabels('')
ax.set_xlabel('')
ax.grid()

ax = fig.add_subplot(324)
ax.plot(Ut, DSu, '-', color='orange')
ax.axis([0,10,0,4])
ax.set_xticklabels('')
ax.set_xlabel('')
ax.grid()

ax = fig.add_subplot(326)
ax.plot(Ut, Qinu/1000, '-m')
ax.axis([0,10,0,70])
ax.set_xlabel('$U_{T}\ (m\ s^{-1})$')
ax.grid()

plt.show()

# RC CLEANUP
plt.rcdefaults()

    

