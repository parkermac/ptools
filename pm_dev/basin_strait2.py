"""
Code to simulate a simple basin and sill system.

This version has TWO BASINS.
"""

import numpy as np
import matplotlib.pyplot as plt
import dev_functions as dfun

testing = True

if testing == False:
    ND = 40 # number of days for integration
else:
    ND = 3*365 # number of days for integration
    
# forcing functions
ff_list = ['linear_ramp_Qr', 'linear_ramp_Ut', 'jumps', 'cycles', 'steady']

# USER: run different cases by changing the number below,
# e.g. ff_list[0] does the first example "linear_ramp_Qr"
force_flag = ff_list[2]


# constants
beta = 7.7e-4
g = 9.8
Cd = 2.6e-3
Socn = 32

# sizes
H = 50 # depth of the sill
B = 3e3 # width of the sill
L = 20e3 # length of the sill
V1 = 50e9 # volume of upper basin box
V2 = 170e9 - V1 # volume of deeper basin box
V3 = (50e9)/3 # volume of upper basin box
V4 = (170e9 - V1)/3 # volume of deeper basin box
# total volume of Puget Sound is about 170 km^3

# derived quantities
c2 = g * beta * H * Socn
# in a more complex system with multiple basins we would
# have Socn be determined dynamically


def fQr(t, force_flag):
    td = t/86400
    if force_flag == 'linear_ramp_Qr':
        Qr = td * 10 + 100
    elif force_flag == 'linear_ramp_Ut':
        Qr = 1000
    elif force_flag == 'jumps':
        Qr = 1000
        if td>=400 and td<500:
            Qr = 2000
        elif td>=500:
            Qr = 500
    elif force_flag == 'cycles':
        ny = np.floor(td/365)
        Qr = 500 + (1000*np.exp(-(td-365/2-ny*365)**2/50**2)
            + 1500*np.exp(-(td-300-ny*365)**2/20**2))
    elif force_flag == 'steady':
        Qr = 1000
    return(Qr)
    
def fUt(t, force_flag):
    td = t/86400
    if force_flag == 'linear_ramp_Qr':
        Ut = 1
    elif force_flag == 'linear_ramp_Ut':
        Ut = .3 + td * (3/700)
    elif force_flag == 'jumps':
        Ut = 1
        if td>=50 and td<250:
            Ut = 1.5
        elif td>=250:
            Ut = .5
    elif force_flag == 'cycles':
        Ut = 1.5 + 0.5*np.cos((2*np.pi/14)*td)*(1 + 0.5*np.cos((4*np.pi/365)*td))
    elif force_flag == 'steady':
        Ut = 1
    return Ut
   
# functions for derived quantities

def fK(Cd, Ut, H):
    a0 = 0.028
    K = a0 * Cd * Ut * H
    Ks = K/2.2
    return K, Ks
    
def fSQin(H, B, c2, K, Ks, L, Qr, Sout_l, Sin_s):
    A = H*B
    H2 = H**2
    L2 = L**2
    a = A*c2*H2/(4*48*Sin_s*K)
    b = 3*H2*a/(20*Ks*A)
    # find dsdx
    dsdx = (-L + np.sqrt(L2 + 8*b*(Sin_s-Sout_l)))/(4*b)
    dq = a * dsdx
    Qin = dq
    Qout = -(Qin + Qr)
    # adjust salinitied to conserve salt flux
    ds = b * dsdx**2
    eps = Qr*(Sin_s-Sout_l-2*ds)/(2*Qin+Qr)
    Sin_l = Sout_l + 2*ds - eps
    Sout_s = Sin_s - 2*ds - eps
    return Sin_l, Sout_s, Qin, Qout, dsdx

# prepare result vectors
    
dt = 86400 # time step (s)
NT = int(ND*dt / 86400)

T = np.linspace(0, NT*dt, NT)
Td = T/86400

S1_a = np.nan * np.ones(NT)
S2_a = np.nan * np.ones(NT)
S3_a = np.nan * np.ones(NT)
S4_a = np.nan * np.ones(NT)
Qin_a = np.nan * np.ones(NT)
Sin_a = np.nan * np.ones(NT)
Qin_ap = np.nan * np.ones(NT) # the landward sill
Sin_ap = np.nan * np.ones(NT) # the landward sill

Qr_a = np.nan * np.ones(NT)
Ut_a = np.nan * np.ones(NT)

# collect forcing time series
for nt in range(NT):
    Qr_a[nt] = fQr(T[nt], force_flag)
    Ut_a[nt] = fUt(T[nt], force_flag)

# draft intial conditions
S1_a[0] = Socn-2
S2_a[0] = Socn-1
S3_a[0] = Socn-3
S4_a[0] = Socn-2

# make an equilibrated initial condition
for nt in range(NT-1):
    Qr = Qr_a[0]
    Ut = Ut_a[0]
    K, Ks = fK(Cd, Ut, H)
    # landward basin
    Sout_l = S3_a[nt]
    Sin_s = S2_a[nt]
    Sin_lp, Sout_sp, Qinp, Qoutp, dsdxp = fSQin(H, B, c2, K, Ks, L, Qr, Sout_l, Sin_s)
    Qin_ap[nt] = Qinp
    Sin_ap[nt] = Sin_lp
    S3_a[nt+1] = S3_a[nt] + dt*( S3_a[nt]*Qoutp/V1 + S4_a[nt]*Qinp/V1)
    S4_a[nt+1] = S4_a[nt] + dt*( Sin_lp*Qinp/V2 - S4_a[nt]*Qinp/V2)
    # seaward basin
    Sout_l = S1_a[nt]
    Sin_s = Socn
    Sin_l, Sout_s, Qin, Qout, dsdx = fSQin(H, B, c2, K, Ks, L, Qr, Sout_l, Sin_s)
    Qin_a[nt] = Qin
    Sin_a[nt] = Sin_l
    Q21 = Qin-Qinp
    if Q21 >= 0:
        F21 = Q21*S2_a[nt]
    else:
        F21 = Q21*S1_a[nt]
    S1_a[nt+1] = S1_a[nt] + dt*( S1_a[nt]*Qout/V1
        - Sout_sp*Qoutp/V1
        + F21/V1 )
    S2_a[nt+1] = S2_a[nt] + dt*( Sin_l*Qin/V2 - S2_a[nt]*Qin/V2)

# reset the initial condition
S1_a[0] = S1_a[-1]
S2_a[0] = S2_a[-1]

# time integration
for nt in range(NT-1):
    Qr = Qr_a[nt+1]
    Ut = Ut_a[nt+1]
    K, Ks = fK(Cd, Ut, H)
    # landward basin
    Sout_l = S3_a[nt]
    Sin_s = S2_a[nt]
    Sin_lp, Sout_sp, Qinp, Qoutp, dsdxp = fSQin(H, B, c2, K, Ks, L, Qr, Sout_l, Sin_s)
    Qin_ap[nt] = Qinp
    Sin_ap[nt] = Sin_lp
    S3_a[nt+1] = S3_a[nt] + dt*( S3_a[nt]*Qoutp/V1 + S4_a[nt]*Qinp/V1)
    S4_a[nt+1] = S4_a[nt] + dt*( Sin_lp*Qinp/V2 - S4_a[nt]*Qinp/V2)
    # seaward basin
    Sout_l = S1_a[nt]
    Sin_s = Socn
    Sin_l, Sout_s, Qin, Qout, dsdx = fSQin(H, B, c2, K, Ks, L, Qr, Sout_l, Sin_s)
    Qin_a[nt] = Qin
    Sin_a[nt] = Sin_l
    Q21 = Qin-Qinp
    if Q21 >= 0:
        F21 = Q21*S2_a[nt]
    else:
        F21 = Q21*S1_a[nt]
    S1_a[nt+1] = S1_a[nt] + dt*( S1_a[nt]*Qout/V1
        - Sout_sp*Qoutp/V1
        + F21/V1 )
    S2_a[nt+1] = S2_a[nt] + dt*( Sin_l*Qin/V2 - S2_a[nt]*Qin/V2)

# finishing up
K, Ks = fK(Cd, Ut, H)
Sout_l = S1_a[-1]
Sin_s = Socn
Sin_l, Sout_s, Qin, Qout, dsdx = fSQin(H, B, c2, K, Ks, L, Qr, Sout_l, Sin_s)
Qin_a[nt+1] = Qin
Sin_a[nt+1] = Sin_l

if False:
    # Checking on salt flux conservation in the strait
    # calculate the net salt flux at both ends of the strait
    # RESULT: it works perfectly now.
    Fl = Qin*Sin_l + Qout*Sout_l
    Fs = Qin*Sin_s + Qout*Sout_s
    print('F landward end = ' + str(Fl))
    print('F seaward end = ' + str(Fs))

# collect forcing time series
for nt in range(NT):
    t = T[nt]
    Qr_a[nt] = fQr(t, force_flag)
    Ut_a[nt] = fUt(t, force_flag)
    
# PLOTTING

def add_line(ax, nd):
    aa = ax.axis()
    ax.plot([nd, nd], aa[-2:], '-k')
    
plt.close()
fig = plt.figure(figsize=(12,8))

ax = fig.add_subplot(311)
ax.plot(Td, S1_a,'-b', Td, S2_a,'-r',
     Td, S3_a, '--b', S4_a, '--r', linewidth=3)
ax.text(.3, .5, '$S_{upper}$', horizontalalignment='left',
    transform=ax.transAxes, color='b', fontsize=20)
ax.text(.3, .3, '$S_{lower}$', horizontalalignment='left',
    transform=ax.transAxes, color='r', fontsize=20)
# ax.text(.1, .5, '$S_{3}$', horizontalalignment='left',
#     transform=ax.transAxes, color='b', fontsize=20)
# ax.text(.1, .3, '$S_{4}$', horizontalalignment='left',
#     transform=ax.transAxes, color='r', fontsize=20)
ax.set_xlim((Td[0],Td[-1]))
ax.set_ylim(32,26)
ax.set_title(force_flag.upper())
ax.grid()

ax = fig.add_subplot(312)
ax.plot(Td, Qin_a/1000, '-m',Td, Qin_ap/1000, '--m', linewidth=3)
ax.text(.1, .8, '$Q_{IN}/1000 (m^3s^{-1})$', horizontalalignment='left',
    transform=ax.transAxes, color='m', fontsize=20)
ax.set_xlim((Td[0],Td[-1]))
aa = ax.axis()
#ax.set_ylim((0, aa[3]))
ax.grid()

ax = fig.add_subplot(313)
ax.plot(Td,Qr_a/1000,'-g', Td,Ut_a,'-c', linewidth=3)
ax.set_xlabel('Time (days)')
ax.set_xlim((Td[0],Td[-1]))
aa = ax.axis()
ax.set_ylim((0,aa[3]))
ax.grid()
ax.text(.1, .8, '$Q_{R}/1000 (m^3s^{-1})$', horizontalalignment='left',
    transform=ax.transAxes, color='g', fontsize=20)
ax.text(.9, .8, '$U_{T} (ms^{-1})$', horizontalalignment='right',
    transform=ax.transAxes, color='c', fontsize=20)

plt.show()
    

