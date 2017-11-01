"""
Code to simulate a simple estuarine channel system.

"""

import numpy as np
import matplotlib.pyplot as plt

testing = True

if testing == True:
    ND = 3 # number of days for integration
else:
    ND = 1000 # number of days for integration
    
# forcing functions
ff_dict = {1:'linear_ramp_Qr',
        2:'linear_ramp_Ut',
        3:'jumps',
        4:'cycles',
        5:'steady',
        6:'Qr_jump'}
# USER: run different cases by changing the number below,
# e.g. ff_list[0] does the first example "linear_ramp_Qr"
force_flag = ff_dict[5]

# constants
beta = 7.7e-4
g = 9.8
Cd = 2.6e-3
Socn = 32

# SIZES
# seaward sill
H = 15 # depth
B = 3e3 # width
L = 200e3 # length

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
    
def dsdq(H, B, K, Ks, dsdx, Socn):
    #print('c2=' + str(c2))
    
    A = H*B
    H2 = H**2
    a = A*c2*H2/(4*48*Socn*K)
    b = 3*H2*a/(20*Ks*A)
    dq = a * dsdx
    ds = b * dsdx**2
    return ds, dq
    
def advance_s(Qr, Ut, dims, sbar, Socn):
    H, B, L, c2, Cd, dt, dx = dims
    NXM = len(sbar)
    NX = NXM+1
    dsdx = np.nan * np.ones(NX)
    dsdx[1:-1] = np.diff(sbar)/dx
    dsdx[0] = dsdx[1]
    dsdx[-1] = dsdx[-2]
    sb = np.nan * np.ones(NX)
    sb[1:-1] = sbar[:-1] + np.diff(sbar)/2
    sb[0] = sb[1]
    sb[-1] = sb[-2] + dsdx[-1]*dx
    K, Ks = fK(Cd, Ut, H)
    ds, dq = dsdq(H, B, K, Ks, dsdx, Socn)
    # qin = dq
    # qout = -(qin + Qr)
    # sin = sb + ds
    # sout = sb - ds
    #sb[-1] = Socn# - ds[-1]
    # fin = qin * sin
    # fout = qout * sout
    # fnet = fin + fout
    fnet = 2*ds*dq - Qr*(sb - ds)
    sbarn = sbar + (dt/(H*B*dx))*np.diff(fnet)
    
    return sbarn
    
# prepare result vectors, arrays
dt = 1800 # time step (s)
NT = int(ND*86400/dt) + 1
T = np.linspace(0, (NT-1)*dt, NT)
Td = T/86400
#
dx = 5e3
NX = int(L/dx)
NXM = NX-1
x = np.linspace(-L, 0, NX)
xm = x[:-1] + np.diff(x)/2
#
Sbar = np.nan * np.ones((NT, NXM))
SBAR = np.nan * np.ones((ND+1, NXM))

# make forcing time series
QR = np.nan * np.ones(NT)
UT = np.nan * np.ones(NT)
for nt in range(NT):
    QR[nt] = fQr(T[nt], force_flag)
    UT[nt] = fUt(T[nt], force_flag)

c2 = g * beta * H * Socn
dims = (H, B, L, c2, Cd, dt, dx)

# draft intial conditions
Sbar[0,:] = Socn * np.exp(xm/20e3)
SBAR[0,:] = Sbar[0,:]

# make an equilibrated initial condition
ii = 1
for nt in range(NT-1):
    sbar = Sbar[nt,:]
    Sbar[nt+1,:] = advance_s(QR[0], UT[0], dims, sbar, Socn)
    if np.mod(Td[nt],1.) == 0:
        SBAR[ii,:] = Sbar[nt+1,:]
        ii+=1
        
    

# PLOTTING
if True:    
    plt.close('all')
    fig = plt.figure(figsize=(12,8))

    ax = fig.add_subplot(111)
    for nd in range(ND):
        ax.plot(xm, SBAR[nd,:])
    ax.set_ylim(-10, 40)

    plt.show()