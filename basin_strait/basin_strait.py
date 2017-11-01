"""
Code to simulate a simple basin and sill system.

This version has TWO BASINS.
"""

import numpy as np
import matplotlib.pyplot as plt

testing = False

if testing == True:
    ND = 40 # number of days for integration
else:
    ND = 3*365 # number of days for integration
    
# forcing functions
ff_dict = {1:'linear_ramp_Qr',
        2:'linear_ramp_Ut',
        3:'jumps',
        4:'cycles',
        5:'steady',
        6:'Qr_jump'}
# USER: run different cases by changing the number below,
# e.g. ff_list[0] does the first example "linear_ramp_Qr"
force_flag = ff_dict[6]

landward_basin = False

# constants
beta = 7.7e-4
g = 9.8
Cd = 2.6e-3
Socn = 32

# SIZES
# seaward sill
H = 50 # depth
B = 3e3 # width
L = 20e3 # length
# seaward basin
V1 = 50e9 # volume of upper basin box
V2 = 170e9 - V1 # volume of deeper basin box
# landward sill
Hp = 50 # depth
Bp = 3e3 # width
Lp = 20e3 # length
# landward basin
V3 = (50e9)/3 # volume of upper basin box
V4 = (170e9 - V1)/3 # volume of deeper basin box
# total volume of Puget Sound is about 170 km^3
dims = (H, B, L, Hp, Bp, Lp, V1, V2, V3, V4)

# forcing functions

def fQr(t, force_flag):
    td = t/86400
    if force_flag == 'linear_ramp_Qr':
        Qr = td * 10 + 100
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
    elif force_flag == 'Qr_jump':
        Qr = 1000
        if td>=100 and td<500:
            Qr = 2000
    else:
        Qr = 1000 # default
    return Qr
    
def fUt(t, force_flag):
    td = t/86400
    if force_flag == 'linear_ramp_Ut':
        Ut = .3 + td * (3/700)
    elif force_flag == 'jumps':
        Ut = 1
        if td>=50 and td<250:
            Ut = 1.5
        elif td>=250:
            Ut = .5
    elif force_flag == 'cycles':
        Ut = 1.5 + 0.5*np.cos((2*np.pi/14)*td)*(1 + 0.5*np.cos((4*np.pi/365)*td))
    else:
        Ut = 1 # default
    return Ut
   
# functions for derived quantities

def fK(Cd, Ut, H):
    a0 = 0.028
    K = a0 * Cd * Ut * H
    Ks = K/2.2
    return K, Ks
    
def fSQin(H, B, K, Ks, L, Qr, Sout_l, Sin_s):
    c2 = g * beta * H * Sin_s
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
    
def advance_s(Qr, Ut, dims, S1, S2, S3, S4, do_check=False, landward_basin=True):
    H, B, L, Hp, Bp, Lp, V1, V2, V3, V4 = dims
    # landward basin
    K, Ks = fK(Cd, Ut, Hp)
    Sout_l = S3
    Sin_s = S2
    Sin_lp, Sout_sp, Qinp, Qoutp, dsdxp = fSQin(Hp, Bp, K, Ks, Lp, Qr, Sout_l, Sin_s)
    S3n = S3 + dt*( S3*Qoutp/V1 + S4*Qinp/V1)
    S4n = S4 + dt*( Sin_lp*Qinp/V2 - S4*Qinp/V2)
    if landward_basin == False:
        S3n = 0
        S4n = 0
        Sout_sp = 0
        Qoutp = Qr
        Qinp = 0
    # seaward basin
    K, Ks = fK(Cd, Ut, H)
    Sout_l = S1
    Sin_s = Socn
    Sin_l, Sout_s, Qin, Qout, dsdx = fSQin(H, B, K, Ks, L, Qr, Sout_l, Sin_s)
    Q21 = Qin-Qinp
    if Q21 >= 0:
        F21 = Q21*S2
    else:
        F21 = Q21*S1
    S1n = S1 + dt*( S1*Qout/V1
        - Sout_sp*Qoutp/V1
        + F21/V1 )
    S2n = S2 + dt*( Sin_l*Qin/V2 - S2*Qin/V2)
    if do_check:
        # Checking on salt flux conservation at both ends
        # of the seaward strait
        # RESULT: it works perfectly now.
        Fl = Qin*Sin_l + Qout*Sout_l
        Fs = Qin*Sin_s + Qout*Sout_s
        print('F landward end = ' + str(Fl))
        print('F seaward end = ' + str(Fs))
    return S1n, S2n, S3n, S4n, Qin, Qinp
    
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
    s1, s2, s3, s4, Qin[nt], Qinp[nt] = advance_s(
            QR[0], UT[0], dims, s1, s2, s3, s4,
            landward_basin=landward_basin)
# reset the initial condition to equilibrated state
S1[0] = s1; S2[0] = s2; S3[0] = s3; S4[0] = s4
#
# actual time integration
for nt in range(NT-1):
    S1[nt+1], S2[nt+1], S3[nt+1], S4[nt+1], Qin[nt], Qinp[nt] = advance_s(
            QR[nt+1], UT[nt+1], dims, S1[nt], S2[nt], S3[nt], S4[nt],
            landward_basin=landward_basin)
#
# finishing up (we already have all the S1-4 final points)
nt = NT-1
junk, junk, junk, junk, Qin[nt], Qinp[nt] = advance_s(
        QR[nt], UT[nt], dims, S1[nt], S2[nt], S3[nt], S4[nt],
        landward_basin=landward_basin)
    
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
    

