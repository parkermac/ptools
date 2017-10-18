"""
Code to simulate a simple basin and sill system.
"""

import numpy as np
import matplotlib.pyplot as plt

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
# total volume of Puget Sound is about 170 km^3

# derived quantities
c2 = g * beta * H * Socn
# in a more complex system with multiple basins we would
# have Socn be determined dynamically

# forcing functions
ff_list = ['linear_ramp_Qr', 'test', 'linear_ramp_Ut', 'jumps', 'cycles']

# USER: run different cases by changing the number below,
# e.g. ff_list[0] does the first example "linear_ramp_Qr"
force_flag = ff_list[3]

def fQr(t, force_flag):
    td = t/86400
    if force_flag in ['linear_ramp_Qr', 'test']:
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
    return(Qr)
    
def fUt(t, force_flag):
    td = t/86400
    if force_flag in ['linear_ramp_Qr', 'test']:
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
        Ut = 1.5 + 0.5*np.cos((2*np.pi/14)*td)*(1 + 0.5*np.cos((2*np.pi/365)*td))
    return Ut
   
# functions for derived quantities

def fK(Cd, Ut, H):
    a0 = 0.028
    K = a0 * Cd * Ut * H
    return K
    
def fQin(H, B, c2, K, L, S1, Socn):
    dsdx = (Socn - (S1)) / L
    ue = (1/48) * c2 * H**2 * dsdx / (Socn * K)
    Qin = 0.5 * H * B * ue
    return Qin
    
def fSin(H, B, c2, K, L, S1, Socn):
    dsdx = (Socn - (S1)) / L
    ue = (1/48) * c2 * H**2 * dsdx / (Socn * K)
    ds = Socn * (3/20) * ue * H**2 * dsdx / (Socn * (K/2.2))
    Sin = S1 + ds
    if Sin >= Socn:
        # giant hack
        Sin = Socn
    return Sin

# prepare result vectors

ND = 2*365 # number of days for integration
#ND = 40 # number of days for integration
dt = 86400 # time step (s)
NT = int(ND*dt / 86400)

T = np.linspace(0, NT*dt, NT)
Td = T/86400
S1 = np.nan * np.ones(NT)
S2 = np.nan * np.ones(NT)
Qin_a = np.nan * np.ones(NT)
Sin_a = np.nan * np.ones(NT)
Qr_a = np.nan * np.ones(NT)
Ut_a = np.nan * np.ones(NT)

# draft intial conditions
S1[0] = Socn - 1
S2[0] = Socn

# make an equilibrated initial condition
for nt in range(NT-1):    
    t = T[0]    
    Qr = fQr(t, force_flag)
    Ut = fUt(t, force_flag)
    K = fK(Cd, Ut, H)
    Qin = fQin(H, B, c2, K, L, S1[nt], Socn)
    Qin_a[nt] = Qin
    Sin = fSin(H, B, c2, K, L, S1[nt], Socn)
    Sin_a[nt] = Sin
    S1[nt+1] = S1[nt] + dt*( - S1[nt]*(Qin + Qr)/V1 + S2[nt]*Qin/V1)    
    S2[nt+1] = S2[nt] + dt*( Sin*Qin/V2 - S2[nt]*Qin/V2)
    
# reset the initial condition
# intial conditions
S1[0] = S1[-1]
S2[0] = S2[-1]

# time integration
for nt in range(NT-1):    
    t = T[nt+1]    
    Qr = fQr(t, force_flag)
    Ut = fUt(t, force_flag)
    K = fK(Cd, Ut, H)
    Qin = fQin(H, B, c2, K, L, S1[nt], Socn)
    Qin_a[nt] = Qin
    Sin = fSin(H, B, c2, K, L, S1[nt], Socn)
    Sin_a[nt] = Sin    
    S1[nt+1] = S1[nt] + dt*( - S1[nt]*(Qin + Qr)/V1 + S2[nt]*Qin/V1)    
    S2[nt+1] = S2[nt] + dt*( Sin*Qin/V2 - S2[nt]*Qin/V2)
# finishing up
Qin = fQin(H, B, c2, K, L, S1[nt+1], Socn)
Qin_a[nt+1] = Qin
Sin = fSin(H, B, c2, K, L, S1[nt+1], Socn)
Sin_a[nt+1] = Sin

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
ax.plot(Td, S1,'-b', Td, S2,'-r', Td, Sin_a, '-c', linewidth=3)
ax.text(.1, .5, '$S_{1}$', horizontalalignment='left',
    transform=ax.transAxes, color='b', fontsize=20)
ax.text(.2, .5, '$S_{2}$', horizontalalignment='left',
    transform=ax.transAxes, color='r', fontsize=20)
ax.text(.3, .5, '$S_{in}$', horizontalalignment='left',
    transform=ax.transAxes, color='c', fontsize=20)
ax.set_xlim((Td[0],Td[-1]))
ax.set_title(force_flag.upper())
ax.grid()

ax = fig.add_subplot(312)
ax.plot(Td,Qin_a/1000,'-g', linewidth=3)
ax.text(.1, .8, '$Q_{IN}/1000 (m^3s^{-1})$', horizontalalignment='left',
    transform=ax.transAxes, color='g', fontsize=20)
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
    

