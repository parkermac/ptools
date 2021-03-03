"""
Explore the Chatwin solution in terms of the efflux-reflux coefficients
for the "L_delta" analyis.

"""

import numpy as np
import matplotlib.pyplot as plt

from importlib import reload
import rflx_fun as rfun
reload(rfun)

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
outdir = os.path.abspath('../../ptools_output/rflx') + '/'
Lfun.make_dir(outdir)

# Naming conventions
#  Layers: s = shallow, d = deep

# define bathymetry and layer thickness
B = 3e3     # width (m)
hs = 20     # thickness of shallow layer (m)
hd = 20     # thickness of deep layer (m)
nx = 1000    # number of steps to initially define a channel

# estuary physical parameters
Socn = 30  # ocean salinity (sbar)
ds = 5 # Sin - Sout at the mouth
Qr = 1e3 # river flow [m3/s]

# get the solution at box edges
Sin, Sout, x, L = rfun.get_Sio_chatwin(Qr, Socn, ds, nx)
DS = Sin - Sout

# get a0 and a1
a0, a1 = rfun.a_calc(Sin, Sout)

# calculate transports at box edges from Knudsen
# sign convention: all Q's positive
Qout = Qr*Sin/(Sin - Sout)
Qin = Qr*Sout/(Sin - Sout)

dx = np.diff(x)
NX = len(x)
xm = x[:-1] + dx/2 # x at box centers
X = x/1000
Xm = xm/1000

# box volumes (at box centers)
dvs = B*hs*dx
dvd = B*hd*dx

# calculate Qup and Qdown at box centers
qdown = a0*Qout[:-1]
qdown[0] = np.nan # avoid divide by zero error
qup = a1*Qin[1:]
ts = dvs/qdown
td = dvd/qup

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


ax1 = fig.add_subplot(311)
# add final model state
ax1.plot(X, Sout, '-', c=cout, lw=lw)
ax1.plot(X, Sin, '-', c=cin, lw=lw)
ax1.grid(True)
ax1.set_ylabel('Salinity')
ax1.text(.1, .8, r'$S_{in}$', c=cin, transform=ax1.transAxes, size=1.5*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax1.text(.2, .8, r'$S_{out}$', c=cout, transform=ax1.transAxes, size=1.5*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax1.text(.82, .25, r'$L_{\Delta}$', c='gray', transform=ax1.transAxes, size=3*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax1.set_xlim(0,X[-1])
ax1.set_xticklabels([])
ax1.set_ylim(0, 1.1*Sin.max())

ax2 = fig.add_subplot(312)
ax2.grid(True)
ax2.set_xticklabels([])
ax2.text(.25, .2, r'$a_{0}$', c=c0, transform=ax2.transAxes, size=1.3*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax2.text(.22, .7, r'$a_{1}$', c=c1, transform=ax2.transAxes, size=1.3*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax2.set_xlim(0,X[-1])
ax2.set_ylim([0,1])

ax3 = fig.add_subplot(313)
ax3.grid(True)
ax3.set_xlabel('X (km)')
ax3.set_ylabel('Time [days]')
ax3.text(.22, .5, r'$T_{shallow}$', c=c0, transform=ax3.transAxes, size=1.3*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax3.text(.25, .05, r'$T_{deep}$', c=c1, transform=ax3.transAxes, size=1.3*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax3.set_xlim(0,X[-1])

# find segments of length L_delta
def find_nearest_ind(array, value):
    # gives the index of the item in array that is closest to value
    idx = (np.abs(array-value)).argmin()
    return idx
# initial segment (closest to mouth)
i0 = NX
i1 = 1
s = Sout[-1]
i0 = find_nearest_ind(Sin, s)
i1 = find_nearest_ind(Sout, s)
while i0 > 0:
    # calculate vertical transports and flushing times over all segments
    # of length L_delta
    ds0 = DS[i0]
    ds1 = DS[i1]
    ds = (ds0 + ds1)/2
    LD = x[i1] - x[i0]
    A0 = (ds1/ds - ds0*ds1/(s*ds))/2
    A1 = (ds0/ds + ds0*ds1/(s*ds))/2
    # print('A0 = %0.1f, A1 = %0.1f' % (A0, A1))
    
    # calculate flushing times:
    Vs = dvs[i0:i1].sum()
    Vd = dvd[i0:i1].sum()
    Qdown = Qout[i0]*A0
    Qup = Qin[i1]*A1
    Ts = Vs/Qdown
    Td = Vd/Qup
    Qe = (Qout[i0] + Qin[i1])/2
        
    ax1.plot([X[i0], X[i1-1]], [s, s], '-k', lw=2*lw, alpha=.3)
    ax2.plot((X[i0] + X[i1])/2, A0, 'o', c=c0)
    ax2.plot((X[i0] + X[i1])/2, A1, 'o', c=c1)
    ax3.plot((X[i0] + X[i1])/2, Ts/86400, 'o', c=c0)
    ax3.plot((X[i0] + X[i1])/2, Td/86400, 'o', c=c1)
    # all remaining segments moving landward
    s = Sout[i0]
    i0 = find_nearest_ind(Sin, s)
    i1 = find_nearest_ind(Sout, s)
    

ax2.set_ylim(bottom=0)
ax3.set_ylim(bottom=0)


# add horizontal flushing time
V = dvs.sum() + dvd.sum()
Th = V/Qout[-1]
ax3.text(.9, .85, r'$T_{flush}=V/Q_{out}=%0.1f\ days$' % (Th/86400), ha='right',
    transform=ax3.transAxes, size=1.3*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))

fig.tight_layout()
fig.savefig(outdir + 'chatwin_a0a1.png')
plt.show()
plt.rcdefaults()
    




