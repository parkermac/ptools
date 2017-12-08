"""
Code to simulate a simple basin and sill system.

This version has TWO BASINS.
"""

import numpy as np
import matplotlib.pyplot as plt
import basin_strait_functions as bsf
from importlib import reload
reload(bsf)

testing = False

if testing:
    ND = 20
    do_check = True
else:
    ND = 1000
    do_check = False

# forcing functions
ff_dict = {1:'linear_ramp_Qr',
        2:'linear_ramp_Ut',
        3:'jumps',
        4:'cycles',
        5:'steady',
        6:'Qr_jump',
        7:'Qr_jump2'}
# USER: run different cases by changing the number below,
# e.g. ff_list[0] does the first example "linear_ramp_Qr"
force_flag = ff_dict[7]

landward_basin = True

dims, params = bsf.get_dims()

# forcing functions

def fQr(t, force_flag):
    td = t/86400
    if force_flag == 'linear_ramp_Qr':
        Qr = td * 10 + 100
        Qr_p = Qr/10
    elif force_flag == 'jumps':
        Qr = 1000
        if td>=400 and td<500:
            Qr = 2000
        elif td>=500:
            Qr = 500
        Qr_p = Qr/10
    elif force_flag == 'cycles':
        ny = np.floor(td/365)
        Qr = 500 + (1000*np.exp(-(td-365/2-ny*365)**2/50**2)
            + 1500*np.exp(-(td-300-ny*365)**2/20**2))
        Qr_p = Qr/10
    elif force_flag == 'Qr_jump':
        Qr = 1000
        if td>=100:# and td<500:
            Qr = 2000
        Qr_p = Qr/10
    elif force_flag == 'Qr_jump2':
        qrp_ratio = False
        Qr = 1000
        Qr_p = 1000
        if td>=100:# and td<500:
            Qr_p = 100
    else:
        Qr = 1000 # default
        Qr_p = Qr/10
    return Qr, Qr_p
    
def fUt(t, force_flag):
    td = t/86400
    if force_flag == 'linear_ramp_Ut':
        Ut = .3 + td * (3/700)
        Ut_p = Ut
    elif force_flag == 'jumps':
        Ut = 1
        if td>=50 and td<250:
            Ut = 1.5
        elif td>=250:
            Ut = .5
        Ut_p = Ut
    elif force_flag == 'cycles':
        Ut = 1.5 + 0.5*np.cos((2*np.pi/14)*td)*(1 + 0.5*np.cos((4*np.pi/365)*td))
        Ut_p = Ut
    else:
        Ut = 1 # default
        Ut_p = Ut
    return Ut, Ut_p
    
def fSo(t, force_flag):
    td = t/86400
    if td >= 400:
        Socn = 22
        Ssfc = 21
    else:
        Socn = 32
        Ssfc = 31
    return Socn, Ssfc
    
# prepare result vectors
    
dt = 86400 # time step (s)
NT = int(ND*dt / 86400)

T = np.linspace(0, NT*dt, NT)
Td = T/86400

S1 = np.nan * np.ones(NT)
S2 = np.nan * np.ones(NT)
S1_p = np.nan * np.ones(NT)
S2_p = np.nan * np.ones(NT)
# landward sill
Qbot = np.nan * np.ones(NT)
# seaward sill
Qbot_p = np.nan * np.ones(NT)

QR = np.nan * np.ones(NT)
UT = np.nan * np.ones(NT)
QR_p = np.nan * np.ones(NT)
UT_p = np.nan * np.ones(NT)
SOCN = np.nan * np.ones(NT)
SSFC = np.nan * np.ones(NT)

# make forcing time series
for nt in range(NT):
    QR[nt], QR_p[nt] = fQr(T[nt], force_flag)
    UT[nt], UT_p[nt] = fUt(T[nt], force_flag)
    SOCN[nt], SSFC[nt] = fSo(T[nt], force_flag)

# draft intial conditions
Socn = SOCN[0]
s1 = Socn-2; s2 = Socn-1; s1_p = Socn-3; s2_p = Socn-2
# make an equilibrated initial condition
rev = False
rev_p = False
if testing == False:
    for nt in range(NT-1):
        s1, s2, s1_p, s2_p, junk, junk, rev, rev_p = bsf.advance_s(
            QR[0], QR_p[0], UT[0], UT_p[0], SOCN[0], SSFC[0],
            dims, params, s1, s2, s1_p, s2_p, dt,
            do_check=do_check, landward_basin=landward_basin, prev_rev=rev, prev_rev_p=rev_p)
# reset the initial condition to equilibrated state
S1[0] = s1; S2[0] = s2; S1_p[0] = s1_p; S2_p[0] = s2_p
#
# actual time integration
for nt in range(NT-1):
        
    S1[nt+1], S2[nt+1], S1_p[nt+1], S2_p[nt+1], Qbot[nt], Qbot_p[nt], rev, rev_p = bsf.advance_s(
        QR[nt+1], QR_p[nt+1], UT[nt+1], UT_p[nt+1], SOCN[nt+1], SSFC[nt+1],
        dims, params,
        S1[nt], S2[nt], S1_p[nt], S2_p[nt], dt,
        do_check=do_check, landward_basin=landward_basin, prev_rev=rev, prev_rev_p=rev_p)
#
# finishing up (we already have all the S1-4 final points)
nt = NT-1
junk, junk, junk, junk, Qbot[nt], Qbot_p[nt], rev, rev_p = bsf.advance_s(
        QR[nt], QR_p[nt], UT[nt], UT_p[nt], SOCN[nt], SSFC[nt],
        dims, params,
        S1[nt], S2[nt], S1_p[nt], S2_p[nt], dt,
        do_check=do_check, landward_basin=landward_basin, prev_rev=rev, prev_rev_p=rev_p)
    
# PLOTTING

    
plt.close()
fig = plt.figure(figsize=(12,8))

ax = fig.add_subplot(311)
ax.plot(Td, S1,'-b', Td, S2,'-r',
     Td, S1_p, '--b', S2_p, '--r', linewidth=3)
ax.text(.05, .85, '$S_{upper}$', horizontalalignment='left',
    transform=ax.transAxes, color='b', fontsize=20)
ax.text(.05, .65, '$S_{lower}$', horizontalalignment='left',
    transform=ax.transAxes, color='r', fontsize=20)
ax.set_xlim((Td[0],Td[-1]))
#ax.set_ylim(32,22)
#ax.set_title(force_flag.upper())
ax.set_title('Seaward Basin/Sill = Solid Line, Landward = Dashed')
ax.grid()

ax = fig.add_subplot(312)
ax.plot(Td, Qbot/1000, '-m',Td, Qbot_p/1000, '--m', linewidth=3)
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
    

