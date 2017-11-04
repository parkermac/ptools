"""
Adjustment time for the basin strait model.
"""

import numpy as np
import matplotlib.pyplot as plt
import basin_strait_functions as bsf
from importlib import reload
reload(bsf)


# forcing functions
ff_dict = {1:'Qr_jump', 2:'Ut_jump'}
force_flag = ff_dict[1]

dims = bsf.get_dims()
H, B, L, Hp, Bp, Lp, V1, V2, V3, V4, beta, g, Cd, Socn = dims

landward_basin = True
ND = 1000 # number of days

# forcing functions

def fQr(t, force_flag):
    td = t/86400
    if force_flag == 'Qr_jump':
        Qr = 1000
        if td>=100:
            Qr = 2000
    else:
        Qr = 1000 # default
    return Qr
    
def fUt(t, force_flag):
    td = t/86400
    if force_flag == 'Ut_jump':
        Ut = 1
        if td>=100:
            Ut = 1.5
    else:
        Ut = 1 # default
    return Ut
    
# prepare result vectors
    
dt = 86400 # time step (s)
NT = int(ND*dt / 86400)
T = np.linspace(0, NT*dt, NT)
Td = T/86400
S1 = np.nan * np.ones(NT)
S2 = np.nan * np.ones(NT)
S3 = np.nan * np.ones(NT)
S4 = np.nan * np.ones(NT)
# landward sill
Qin = np.nan * np.ones(NT)
Sin = np.nan * np.ones(NT)
# seaward sill
Qinp = np.nan * np.ones(NT)
Sinp = np.nan * np.ones(NT)
# forcing
QR = np.nan * np.ones(NT)
UT = np.nan * np.ones(NT)
# make forcing time series
for nt in range(NT):
    QR[nt] = fQr(T[nt], force_flag)
    UT[nt] = fUt(T[nt], force_flag)

# draft intial conditions
s1 = Socn-2; s2 = Socn-1; s3 = Socn-3; s4 = Socn-2
# make an equilibrated initial condition
for nt in range(NT-1):
    s1, s2, s3, s4, Qin[nt], Qinp[nt] = bsf.advance_s(
            QR[0], UT[0], dims, s1, s2, s3, s4, dt,
            landward_basin=landward_basin)
# reset the initial condition to equilibrated state
S1[0] = s1; S2[0] = s2; S3[0] = s3; S4[0] = s4
#
# actual time integration
for nt in range(NT-1):
    S1[nt+1], S2[nt+1], S3[nt+1], S4[nt+1], Qin[nt], Qinp[nt] = bsf.advance_s(
            QR[nt+1], UT[nt+1], dims, S1[nt], S2[nt], S3[nt], S4[nt], dt,
            landward_basin=landward_basin)
#
# finishing up (we already have all the S1-4 final points)
nt = NT-1
junk, junk, junk, junk, Qin[nt], Qinp[nt] = bsf.advance_s(
        QR[nt], UT[nt], dims, S1[nt], S2[nt], S3[nt], S4[nt], dt,
        landward_basin=landward_basin)

# calculate the Qin adjustment time
def calc_tadj(ser, Td):
    v0 = ser[0]
    v1 = ser[-1]
    dv = v1-v0
    vv = (ser - v0)/dv
    value = (1 - 1/np.e)
    idx_q = (np.abs(vv-value)).argmin()
    idx_t0 = (np.abs(Td-100)).argmin()
    t0 = Td[idx_t0]
    t1 = Td[idx_q]
    tadj = t1-t0
    print('Adjustment Time = %0.2f days' % (tadj))

for vv in [S1, S2, Qin]:
    calc_tadj(vv, Td)
    

# PLOTTING

def add_line(ax, nd):
    aa = ax.axis()
    ax.plot([nd, nd], aa[-2:], '-k')
    
#plt.close()
fig = plt.figure(figsize=(12,8))

ax = fig.add_subplot(311)
ax.plot(Td, S1,'-b', Td, S2,'-r',
     Td, S3, '--b', S4, '--r', linewidth=3)
ax.text(.05, .85, '$S_{upper}$', horizontalalignment='left',
    transform=ax.transAxes, color='b', fontsize=20)
ax.text(.05, .65, '$S_{lower}$', horizontalalignment='left',
    transform=ax.transAxes, color='r', fontsize=20)
ax.set_xlim((Td[0],Td[-1]))
ax.set_ylim(32,22)
#ax.set_title(force_flag.upper())
ax.set_title('Seaward Basin/Sill = Solid Line, Landward = Dashed')
ax.grid()

ax = fig.add_subplot(312)
ax.plot(Td, Qin/1000, '-m',Td, Qinp/1000, '--m', linewidth=3)
ax.text(.95, .8, '$Q_{IN}/1000 (m^3s^{-1})$', horizontalalignment='right',
    transform=ax.transAxes, color='m', fontsize=20)
ax.set_xlim((Td[0],Td[-1]))
aa = ax.axis()
#ax.set_ylim((0, aa[3]))
ax.grid()

ax = fig.add_subplot(313)
ax.plot(Td,QR/1000,'-g', linewidth=6)
ax.plot(Td,UT,'-k', linewidth=2)
ax.set_xlabel('Time (days)')
ax.set_xlim((Td[0],Td[-1]))
aa = ax.axis()
ax.set_ylim((0,aa[3]))
ax.grid()
ax.text(.05, .1, '$Q_{R}/1000 (m^3s^{-1})$', horizontalalignment='left',
    transform=ax.transAxes, color='g', fontsize=20)
ax.text(.95, .8, '$U_{T} (ms^{-1})$', horizontalalignment='right',
    transform=ax.transAxes, color='k', fontsize=20)

plt.show()
    

