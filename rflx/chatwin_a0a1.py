"""
Explore the Chatwin solution in terms of the efflux-reflux coefficients
for the "L_delta" analyis.

"""

import numpy as np
import matplotlib.pyplot as plt

from importlib import reload
import rflx_fun as rfun
reload(rfun)

# Naming conventions
#  Layers: s = shallow, d = deep

# define bathymetry and layer thickness
B = 3e3     # width (m)
hs = 10     # thickness of shallow layer (m)
hd = 20     # thickness of deep layer (m)
nx = 1000    # number of steps to initially define a channel

# estuary physical parameters
Socn = 30  # ocean salinity (sbar)
ds = 5 # Sin - Sout at the mouth
Qr = 1e3 # river flow [m3/s]

# get the solution at box edges
Sin, Sout, x, L = rfun.get_Sio_chatwin(Qr, Socn, ds, nx)
DS = Sin - Sout
X = x/1000

# calculate transports at box edges from Knudsen
# sign convention: all Q's positive
Qout = Qr*Sin/(Sin - Sout)
Qin = Qr*Sout/(Sin - Sout)

dx = np.diff(x)
NX = len(x)
xm = x[:-1] + dx/2 # x at box centers

# box volumes (at box centers)
dvs = B*hs*dx
dvd = B*hd*dx

# find a segment of length LD
def find_nearest_ind(array, value):
    # gives the index of the item in array that is closest to value
    idx = (np.abs(array-value)).argmin()
    return idx

i0 = NX
i1 = 1
s = Sout[-1]
i0 = find_nearest_ind(Sin, s)
i1 = find_nearest_ind(Sout, s)

# PLOTTING
plt.close('all')

fs = 16
lw = 3
plt.rc('font', size=fs)
fig = plt.figure(figsize=(12,8))
cin = 'r'
cout = 'b'
c0 = 'g'
c1 = 'darkorange'

ax1 = fig.add_subplot(211)
# add final model state
ax1.plot(X, Sout, '-', c=cout, lw=lw)
ax1.plot(X, Sin, '-', c=cin, lw=lw)
ax1.grid(True)
ax1.set_ylabel('Salinity')
ax1.text(.1, .8, r'$S_{in}$', c=cin, transform=ax1.transAxes, size=1.5*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax1.text(.2, .8, r'$S_{out}$', c=cout, transform=ax1.transAxes, size=1.5*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax1.text(.8, .3, r'$L_{\Delta}$', c='gray', transform=ax1.transAxes, size=3*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax1.set_xlim(0,X[-1])
ax1.set_ylim(0, 1.1*Sin.max())

ax2 = fig.add_subplot(212)
# ax2.plot(X, a0, '-', c=c0, lw=lw)
# ax2.plot(X, a1, '-', c=c1, lw=lw)
ax2.grid(True)
ax2.set_xlabel('X (km)')
ax2.set_ylabel('Time [days]')
ax2.text(.1, .85, r'$T^{reflux}_{flush}$', c=c0, transform=ax2.transAxes, size=1.3*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax2.text(.1, .6, r'$T^{efflux}_{flush}$', c=c1, transform=ax2.transAxes, size=1.3*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax2.set_xlim(0,X[-1])

while i0 > 0:
    
    ds0 = DS[i0]
    ds1 = DS[i1]
    ds = (ds0 + ds1)/2

    LD = x[i1] - x[i0]

    a0 = (ds1/ds - ds0*ds1/(s*ds))/2
    a1 = (ds0/ds + ds0*ds1/(s*ds))/2

    # calculate flushing times:
    Vs = dvs[i0:i1].sum()
    Vd = dvd[i0:i1].sum()
    Qdown = Qout[i0]*a0
    Qup = Qin[i1]*a1
    Ts = Vs/Qdown
    Td = Vd/Qup
    Qe = (Qout[i0] + Qin[i1])/2

    if False:
        print('\ns = %0.1f' % (s))
        print('Shallow vertical flushing time = %0.1f days' % (Ts/86400))
        print('Deep vertical flushing time = %0.1f days' % (Td/86400))
    
    ax1.plot([X[i0], X[i1-1]], [s, s], '-k', lw=2*lw, alpha=.3)
    ax2.plot((X[i0] + X[i1])/2, Ts/86400, 'o', c=c0)
    ax2.plot((X[i0] + X[i1])/2, Td/86400, 'o', c=c1)
    
    s = Sout[i0]
    i0 = find_nearest_ind(Sin, s)
    i1 = find_nearest_ind(Sout, s)

ax2.set_ylim(bottom=0)

# add horizontal flushing time
V = dvs.sum() + dvd.sum()
Th = V/Qout[-1]
ax2.text(.9, .85, r'$T_{flush}=V/Q_{out}=%0.1f\ days$' % (Th/86400), ha='right',
    transform=ax2.transAxes, size=1.3*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))

plt.show()
plt.rcdefaults()
    




