"""
Solve the efflux-reflux system with a 1-segment Chatwin solution.

This is focused on testing how well the numerical scheme (the flux engine)
reproduced the original salinities.

Using new notation.

"""

import numpy as np
import matplotlib.pyplot as plt

from importlib import reload
import new_rflx_fun as rfun
reload(rfun)

# Naming conventions
#  Layers: s = shallow, d = deep

# define bathymetry and layer thickness
B = 3e3     # width (m)
hs = 10     # thickness of shallow layer (m)
hd = 20     # thickness of deep layer (m)
nx = 100    # number of steps to initially define a channel

# estuary physical parameters
Socn = 30  # ocean salinity
ds = 5 # Sin - Sout at the mouth
Qr = 1e3 # river flow [m3/s]

# get the solution at box edges
Sin, Sout, x, L = rfun.get_Sio_chatwin(Qr, Socn, ds, nx)

# calculate transports at box edges from Knudsen
Qout = -Qr*Sin/(Sin - Sout)
Qin = Qr*Sout/(Sin - Sout)

# run specifications
ndays = 200
dx = np.diff(x)
NX = len(x)
xm = x[:-1] + dx/2 # x at box centers

# box volumes (at box centers)
dvs = B*hs*dx
dvd = B*hd*dx

# calculate Efflux-Reflux fractions (a11 and a00) for each box
# using from,to notation.
a00, a11 = rfun.a_calc(Sin, Sout)

dt, NT, NS = rfun.get_time_step(dx, Qout, B, hs, ndays)

# pack some parameters
info_tup = (NS, NX, NT, dt, dvs, dvd, Qout, Qin, a00, a11)

# intial condition vectors
csp = np.zeros(NX-1) # csp = "concentration shallow previous"
cdp = np.zeros(NX-1) # cdp = "concentration deep previous"

exp = 'ocean_one'
plot_s = False
if exp == 'salt':
    # Reproduce salinity state
    csa, cda = rfun.c_calc(csp, cdp, info_tup, ocn=Sin[-1])
    plot_s = True
elif exp == 'all_one':
    # Conserve a constant tracer to test volume conservation
    csa, cda = rfun.c_calc(csp, cdp, info_tup, riv=1, ocn=1)
elif exp == 'river_one':
    # Constant river source to test tracer conservation
    csa, cda = rfun.c_calc(csp, cdp, info_tup, riv=1, ocn=0)
    print('Tracer flux out at mouth %0.3f [c m3/s]' % (-Qout[-1]*csa[-1,-1]))
    print('Tracer flux in at head %0.3f [c m3/s]' % (1*Qr))
elif exp == 'ocean_one':
    # Constant river source to test tracer conservation
    csa, cda = rfun.c_calc(csp, cdp, info_tup, riv=0, ocn=1)
    print('Tracer flux out at mouth %0.3f [c m3/s]' % (-Qout[-1]*csa[-1,-1]))
    print('Tracer flux in at mouth %0.3f [c m3/s]' % (1*Qin[-1]))
else:
    print('exp = %s not supported' % (exp))
    sys.exit()

# PLOTTING
plt.close('all')

X = x/1e3
Xm = xm/1e3

fig = plt.figure(figsize=(17,10))

ax = fig.add_subplot(211)
# add final model state
ax.plot(Xm, csa[-1,:], '-b')
ax.plot(Xm, cda[-1,:], '-r')
if plot_s:
    ax.plot(X, Sin, ':r', X, Sout, ':b', linewidth=2)
    ax.text(.1, .8, 'Dotted = Target Solution', transform=ax.transAxes)
ax.grid(True)
ax.set_xlabel('X (km)')
ax.set_ylabel('Sin, Sout')
ax.text(.1, .9, 'Red = In, Blue = Out', transform=ax.transAxes)
ax.text(.1, .7, 'Solid = End Numerical Solution', transform=ax.transAxes)
ax.set_xlim(0,X[-1])
ax.set_ylim(0, 1.1*cda[-1,:].max())

ax = fig.add_subplot(212)
ax.plot(Xm, a11, '-r',Xm, a00, '-b', linewidth=2)
ax.grid(True)
ax.set_xlabel('X (km)')
ax.set_ylabel('a11, a00')
ax.text(.1, .9, 'Red = a11, Blue = a00', transform=ax.transAxes)
ax.set_xlim(0,X[-1])
ax.set_ylim(bottom=0)

plt.show()



