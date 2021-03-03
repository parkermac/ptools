"""
Solve the efflux-reflux system with two kinds of estuaries.

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
nx = 100    # number of steps to initially define a channel
# estuary physical parameters
Socn = 30  # ocean salinity
ds = 5 # Sin - Sout at the mouth
Qr = 1e3 # river flow [m3/s]

plt.close('all')

# **** Himix estuary ****

# get the solution at box edges
Sin, Sout, x, L = rfun.get_Sio_chatwin(Qr, Socn, ds, nx)
# calculate transports at box edges from Knudsen
# sign convention: all Q's positive
Qout = Qr*Sin/(Sin - Sout)
Qin = Qr*Sout/(Sin - Sout)
# run specifications
ndays = 200
dx = np.diff(x)
NX = len(x)
xm = x[:-1] + dx/2 # x at box centers
X = x/1e3
Xm = xm/1e3
# box volumes (at box centers)
dvs = B*hs*dx
dvd = B*hd*dx
# calculate Efflux-Reflux fractions (a11 and a00) for each box
a0, a1 = rfun.a_calc(Sin, Sout)
# get time step
dt, NT, NS = rfun.get_time_step(dx, Qout, B, hs, ndays)
# pack some parameters
info_tup = (NS, NX, NT, dt, dvs, dvd, Qout, Qin, a0, a1)
# intial condition vectors
csp = np.zeros(NX-1) # csp = "concentration shallow previous"
cdp = np.zeros(NX-1) # cdp = "concentration deep previous"
# experiments
csa_salt, cda_salt = rfun.c_calc(csp, cdp, info_tup, riv=0, ocn=Sin[-1])
csa_one, cda_one = rfun.c_calc(csp, cdp, info_tup, riv=1, ocn=1)
csa_oxy, cda_oxy = rfun.c_calc(csp, cdp, info_tup, riv=1, ocn=1, Ts=10, Td=10, Cs=1, Cd=0)
# PLOTTING
fs = 16
lw = 3
plt.rc('font', size=fs)
fig = plt.figure(figsize=(12,10))
cin = 'r'
cout = 'b'
c0 = 'g'
c1 = 'darkorange'
# salt
ax = fig.add_subplot(221)
ax.plot(Xm, csa_salt[-1,:], '-', c=cout, lw=lw)
ax.plot(Xm, cda_salt[-1,:], '-', c=cin, lw=lw)
ax.grid(True)
ax.set_xlim(0,X[-1])
ax.set_ylabel('Salinity')
ax.text(.1, .8, '$S_{in}$', c=cin, transform=ax.transAxes, size=1.5*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax.text(.2, .8, '$S_{out}$', c=cout, transform=ax.transAxes, size=1.5*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax.set_ylim(0, 1.1*Sin[-1])
ax.set_xlabel('X (km)')
ax.set_title('High Mixing Estuary')
# tracer
ax = fig.add_subplot(223)
ax.plot(Xm, csa_one[-1,:], '-', c=c0, lw=lw)
ax.plot(Xm, cda_one[-1,:], '-', c=c1, lw=lw)
ax.plot(Xm, csa_oxy[-1,:], '--', c=c0, lw=lw)
ax.plot(Xm, cda_oxy[-1,:], '--', c=c1, lw=lw)
ax.grid(True)
ax.set_xlim(0,X[-1])
ax.set_ylabel('Concentration')
ax.text(.1, .8, '$C_{in}$', c=c1, transform=ax.transAxes, size=1.5*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax.text(.2, .8, '$C_{out}$', c=c0, transform=ax.transAxes, size=1.5*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax.set_ylim(0, 1.1)
ax.set_xlabel('X (km)')

# **** Lowmix estuary ****
ds = (Sin-Sout).mean()
# get the solution at box edges
Sin, Sout, x, L = rfun.get_Sio_fjord(Socn, ds, nx, L=200e3)
# calculate transports at box edges from Knudsen
# sign convention: all Q's positive
Qout = Qr*Sin/(Sin - Sout)
Qin = Qr*Sout/(Sin - Sout)
Qout[0] = Qr
Qin[0] = 0
# run specifications
ndays = 200
dx = np.diff(x)
NX = len(x)
xm = x[:-1] + dx/2 # x at box centers
X = x/1e3
Xm = xm/1e3
# box volumes (at box centers)
dvs = B*hs*dx
dvd = B*hd*dx
# calculate Efflux-Reflux fractions (a11 and a00) for each box
a0, a1 = rfun.a_calc(Sin, Sout)
a0[0] = 0
a1[0] = 1
# get time step
dt, NT, NS = rfun.get_time_step(dx, Qout, B, hs, ndays*3)
# pack some parameters
info_tup = (NS, NX, NT, dt, dvs, dvd*5, Qout, Qin, a0, a1)
# intial condition vectors
csp = np.zeros(NX-1) # csp = "concentration shallow previous"
cdp = np.zeros(NX-1) # cdp = "concentration deep previous"
# experiments
csa_salt, cda_salt = rfun.c_calc(csp, cdp, info_tup, riv=0, ocn=Sin[-1])
csa_one, cda_one = rfun.c_calc(csp, cdp, info_tup, riv=1, ocn=1)
csa_oxy, cda_oxy = rfun.c_calc(csp, cdp, info_tup, riv=1, ocn=1, Ts=10, Td=10, Cs=1, Cd=0)
# PLOTTING
# salt
ax = fig.add_subplot(222)
ax.plot(Xm, csa_salt[-1,:], '-', c=cout, lw=lw)
ax.plot(Xm, cda_salt[-1,:], '-', c=cin, lw=lw)
ax.grid(True)
ax.set_xlim(0,X[-1])
ax.set_ylabel('Salinity')
ax.text(.1, .8, '$S_{in}$', c=cin, transform=ax.transAxes, size=1.5*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax.text(.2, .8, '$S_{out}$', c=cout, transform=ax.transAxes, size=1.5*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax.set_ylim(0, 1.1*Sin[-1])
ax.set_xlabel('X (km)')
ax.set_title('Low Mixing Estuary')
# tracer
ax = fig.add_subplot(224)
ax.plot(Xm, csa_one[-1,:], '-', c=c0, lw=lw)
ax.plot(Xm, cda_one[-1,:], '-', c=c1, lw=lw)
ax.plot(Xm, csa_oxy[-1,:], '--', c=c0, lw=lw)
ax.plot(Xm, cda_oxy[-1,:], '--', c=c1, lw=lw)
ax.grid(True)
ax.set_xlim(0,X[-1])
ax.set_ylabel('Concentration')
ax.text(.1, .8, '$C_{in}$', c=c1, transform=ax.transAxes, size=1.5*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax.text(.2, .8, '$C_{out}$', c=c0, transform=ax.transAxes, size=1.5*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax.set_ylim(0, 1.1)
ax.set_xlabel('X (km)')

fig.tight_layout()
fig.savefig(outdir + 'lowmix_himix.png')
plt.show()

plt.rcdefaults()




