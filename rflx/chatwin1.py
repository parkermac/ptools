"""
Solve the efflux-reflux system with a 1-segment Chatwin solution.
"""

import numpy as np
import matplotlib.pyplot as plt

from importlib import reload
import rflx_fun as rfun
reload(rfun)

# Naming conventions
#  Layers: s = shallow, d = deep
#  Channels: o_ = outer, i_ = inner

B = 3e3     # width (m)
hs = 10     # thickness of shallow layer (m)
hd = 10     # thickness of deep layer (m)
nx = 200    # number of steps to initially define a channel

# parameters
Socn = 30 
ds = 5
Qr = 1e3

# first get the outer solution
Sin, Sout, x, L = rfun.get_Sio_chatwin(Qr, Socn, ds, nx)

# find outer transports
Qin = Qr*Sout/(Sin - Sout)
Qout = -Qr*Sin/(Sin - Sout)

# time dependent run specifications
ndays = 100
dx = np.diff(x)
NX = len(x)
xm = x[:-1] + dx/2

# cell volumes (at bin centers)
dvs = B*hs*dx
dvd = B*hd*dx
# total volume
V = np.sum(dvs + dvd)

# calculate Efflux-Reflux fractions (a21 and a34), at xm
o_Qr = 0
I = 0
q1, q2, q3, q4, a21, a34 = rfun.qa_calc(Qin, Qout, Sin, Sout, I)

# calculate the time step dynamically using some factor (0.9) times
# the shortest time for water to traverse a grid box in the
# along_channel direction (experiments showed that it was unstable
# for cfl_factor > 1.1)
cfl_factor = 0.9
dt = cfl_factor * np.min(dx/(-Qout[1:]/(B*hs))) # time step (s)
NT = int(ndays*86400/dt) # total number of time steps
NS = 10 # number of saves

# pack some parameters
info_tup = (NS, NX, NT, dt, dvs, dvd, q1, q2, q3, q4, a21, a34, o_Qr, I)

# PLOTTING
plt.close('all')

X = x/1e3
Xm = xm/1e3

fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(11,7))

ax = axes[0,0]
ax.plot(X, Sin, ':r', X, Sout, ':b', linewidth=2, alpha=.5)
ax.grid(True)

ax = axes[1,0]
ax.plot(X, Qin/1e3, ':r', X, -Qout/1e3, ':b', linewidth=2, alpha=.5)
ax.grid(True)
    
# reproduce salinity state
# intial condition vectors
csp = np.zeros(NX-1) # csp = "concentration shallow previous"
cdp = np.zeros(NX-1) # cdp = "concentration deep previous"
csa, cda, f_tup = rfun.c_calc(csp, cdp, info_tup, ocn=Sin[-1])

ax = axes[0,0]
ax.plot(Xm, csa[-1,:], '-b')
ax.plot(Xm, cda[-1,:], '-r')

# tracer in river
csa, cda, f_tup = rfun.c_calc(csp, cdp, info_tup, riv=1)

ax = axes[2,0]
ax.plot(Xm, csa[-1,:], '-b')
ax.plot(Xm, cda[-1,:], '-r')

T, c_tot, f1_vec, f2_vec, f3_vec, f4_vec, friv2_vec = f_tup
f_net = f1_vec + f2_vec + f3_vec + f4_vec + friv2_vec
c_tot_from_f = (c_tot[0] +
    np.concatenate(([0.], np.cumsum(f_net*dt)[:-1])))

ax = axes[0,1]
ax.plot(T, c_tot/V, '-c', linewidth=4)
ax.plot(T, c_tot_from_f/V, '-k')

ax = axes[1,1]
ax.plot(T, f4_vec, '-k', label='River Source')
ax.plot(T, -f2_vec, '--k', label='Actual Mouth Sink')
ax.legend()

plt.show()



