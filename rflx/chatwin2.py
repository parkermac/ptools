"""
Solve the efflux-reflux system with a 2-segment Chatwin solution
in which a second river enters midway along the channel.
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
nx = 500    # number of steps to initially define a channel

# outer
o_Socn = 30 
o_ds = 5
o_Qr = 2e3

# inner
i_Qr = 1e3

# first get the outer solution
o_Sin, o_Sout, o_x, o_L = rfun.get_Sio_chatwin(o_Qr + i_Qr, o_Socn, o_ds, nx)

# target x-position at which to match solutions
xio = o_L/3
# and find its index in the outer solution
io = np.argmax(o_x>xio)
# o_x position where solutions match
xio = o_x[io]

# trim outer solution
o_Sin = o_Sin[io:]
o_Sout = o_Sout[io:]
o_x = o_x[io:]

# find outer transports
o_Qin = (o_Qr + i_Qr)*o_Sout/(o_Sin - o_Sout)
o_Qout = -(o_Qr + i_Qr)*o_Sin/(o_Sin - o_Sout)

# then get the inner solution
i_Socn = (o_Sin[0] + o_Sout[0])/2
i_ds = o_Sin[0] - o_Sout[0]
i_Sin, i_Sout, i_x, i_L = rfun.get_Sio_chatwin(i_Qr, i_Socn, i_ds, nx)
# find inner transports
i_Qin = i_Qr*i_Sout/(i_Sin - i_Sout)
i_Qout = -i_Qr*i_Sin/(i_Sin - i_Sout)

# adjust the first x of outer solution to match last x of inner solution
# plus one extra space inbetween
o_x = o_x + i_x[-1] - xio + (i_x[1] - i_x[0])

# concatenate to make a single vector solution
Qin = np.concatenate((i_Qin, o_Qin))
Qout = np.concatenate((i_Qout, o_Qout))
Sin = np.concatenate((i_Sin, o_Sin))
Sout = np.concatenate((i_Sout, o_Sout))
x = np.concatenate((i_x, o_x))

# note the length of the inner solution
# which we use for the index "I" of the river source
nix = len(i_x)
I = nix-1 # location of the junction

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
q1, q2, q3, q4, a21, a34 = rfun.qa_calc(Qin, Qout, Sin, Sout, I)

# calculate the time step dynamically using some factor (0.9) times
# the shortest time for water to traverse a grid box in the
# along_channel direction (experiments showed that it was unstable
# for cfl_factor > 1.1)
cfl_factor = 0.9
dt = cfl_factor * np.min(dx/(-Qout[1:]/(B*hs))) # time step (s)
NT = int(ndays*86400/dt) # total number of time steps
NS = 10 # number of saves

# pack some parameter
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

# tracer in second river
csa, cda, f_tup = rfun.c_calc(csp, cdp, info_tup, riv2=1)

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
ax.plot(T, friv2_vec, '-k', label='Second River Source')
ax.plot(T, -f2_vec, '--k', label='Mouth Sink')
ax.legend()

plt.show()



