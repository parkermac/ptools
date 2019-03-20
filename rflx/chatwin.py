"""
Plot the Chatwin TEF and M solutions.
"""

import numpy as np
import matplotlib.pyplot as plt

# Naming convention: s = shallow, d = deep for the two layers

# estuary physical specification
Socn = 30   # ocean salinity (approximate)
ds = 5      # delta_s at mouth
L=50e3      # length of salt intrusion (m)
B = 3e3     # width (m)
hs = 10     # thickness of upper layer (m)
hd = 10     # thickness of lower layer (m)
Qr = 1e3    # river flow (m3/s)
# eventually hs, hd, B, dx and even Qr should be vectors, although then
# the Chatwin solution would not be valid, and one would instead take
# the Sin and Sout from a numical model analyzed using TEF.

# run specifications
NX = 500
ndays = 200
dx = L/NX
exp='' # experiment name

# cell volumes (at bin centers)
dvs = B*hs*dx * np.ones(NX-1)
dvd = B*hd*dx * np.ones(NX-1)
# total volume
V = np.sum(dvs + dvd)

# sinking parameters
wsink = 10 # particle sinking speed (m/day)
sfact = (wsink/86400)*(B*dx) # sinking factor (m3/s)

# note (dt/dv)*B*dx = dt/h; does that make sense?

if True:
    # Chatwin solution
    #
    # convenience factors that define the solution
    a = Socn/(L**1.5)
    alpha = ds/L
    #
    # define x axes (at bin edges)
    # the lower limit of x ensures that Qin[0] = 0
    x = np.linspace((alpha/(2*a))**2,L,NX)
    #
    # Chatwin solution in terms of TEF quantities
    Sin = a*x**1.5 + alpha*x/2
    Sout = a*x**1.5 - alpha*x/2
    if False:
        # analytical form
        Qin = Qr*(a/alpha)*x**0.5 - Qr/2
        Qout = -(Qr*(a/alpha)*x**0.5 + Qr/2)
    else:
        # using Knudsen
        DS = Sin - Sout
        Qin = Qr*Sout/DS
        Qout = -Qr*Sin/DS
    #
    # Mixing: Volume integral landward of a section of the
    # rate of destruction of variance
    M = Qr*(a**2*x**3 - alpha**2*x**2/4)
    # presumably this is equal to Qr*Sin*Sout
    # It would be interesting to also plot dM/dx * dx, the local mixing
else:
    # A very simple solution, just as a proof of concept.
    # The idea is to specify the salinity fields, and then
    # use Knudsen to do the rest.
    x = np.linspace(0,L,NX)
    Sin = np.linspace(ds/10,Socn,NX)
    Sout = np.linspace(0,Socn-ds,NX)
    DS = Sin - Sout
    Qin = Qr*Sout/DS
    Qout = -Qr*Sin/DS
    
# some other calculations from the target salinity field
socn = Sin[-1] # actual ocean salinity
# these are on the midpoints
sin = Sin[:-1] + np.diff(Sin)/2
sout = Sout[:-1] + np.diff(Sout)/2
# freshwater fraction
fws_frac = (socn - sout)/socn
fwd_frac = (socn - sin)/socn
fwv = np.sum(fws_frac*dvs + fwd_frac*dvd)
tres_fw = fwv/(Qr*86400)

# midpoint positions
xm = x[:-1] + np.diff(x)/2

# calculate Efflux-Reflux fractions (a21 and a34), at xm
q1 = Qin[1:]
q3 = Qin[:-1]
q2 = -Qout[1:]
q4 = -Qout[:-1]
s1 = Sin[1:]
s3 = Sin[:-1]
s2 = Sout[1:]
s4 = Sout[:-1]
# notation follows Cokelet and Stewart (1985)
a21 = (q2/q1)*(s2-s4)/(s1-s4) # efflux - up
a34 = (q3/q4)*(s1-s3)/(s1-s4) # reflux - down

# calculate the time step dynamically using some factor (0.9) times
# the shortest time for water to traverse a grid box in the
# along_channel direction (experimets showed that it was unstable for > 1.1)
dt = 0.9 * np.min(dx/(-Qout/(B*hs))) # time step (s)
bs = dt/dvs # dt/dvol shallow
bd = dt/dvd # dt/dvol deep
NT = int(ndays*86400/dt)
NS = 10 # number of saves
NTs = int(NT/NS) # save interval, in timesteps


def c_calc(riv=0, ocn=0, sfact=0):
    # time integration of a boundary value problem
    # NOTE: all variables defined in the calling program are
    # implicitly GLOBAL and so can be used in this function.
    #
    # intial condition vectors
    if exp == 'salt_ic':
        # start with the tracer equal to the salinity field
        csp = sout
        cdp = sin
    else:
        # start with all zeros
        csp = np.zeros(NX-1)
        cdp = np.zeros(NX-1)
    # arrays to save in
    csa = np.nan * np.ones((NS,NX-1))
    cda = np.nan * np.ones((NS,NX-1))
    # time series to save in
    c_tot = np.nan * np.ones(NT) # net amount of tracer
    t_vec = np.nan * np.ones(NT) # time (s)
    f_net = np.nan * np.ones(NT) # net influx of tracer
    #
    tta = 0
    for tt in range(NT):
        # boundary conditions
        c4 = np.concatenate(([riv], csp[:-1])) # river
        c1 = np.concatenate((cdp[1:], [ocn])) # ocean
        # save time series entries
        t_vec[tt] = tt*dt
        c_tot[tt] = np.sum(csp*dvs + cdp*dvd)
        f_net[tt] = q1[-1]*c1[-1] - q2[-1]*csp[-1] + q4[0]*c4[0] - q3[0]*cdp[0]
        # effect of net fluxes
        cs = csp + bs*(q4*c4 + a21*q1*cdp - a34*q4*csp - q2*csp - sfact*csp)
        cd = cdp + bd*(q1*c1 - a21*q1*cdp + a34*q4*csp - q3*cdp + sfact*csp)
        csp = cs.copy()
        cdp = cd.copy()
        if np.mod(tt, NTs) == 0:
            csa[tta,:] = cs
            cda[tta,:] = cd
            tta += 1
    # work on time series
    td_vec = t_vec/86400
    c_tot_from_f = (c_tot[0] +
        np.concatenate(([0.], np.cumsum(f_net*dt)[:-1])))
    
    return csa, cda, td_vec, c_tot, c_tot_from_f
    
def age_calc(riv=0, ocn=0):
    # returns age of water entering from river or ocean
    #
    # intial condition vectors
    # start with all zeros
    csp = np.zeros(NX-1) # normal tracer
    cdp = np.zeros(NX-1)
    ccsp = np.zeros(NX-1) # aging tracer
    ccdp = np.zeros(NX-1)
    #
    tta = 0
    for tt in range(NT):
        # boundary conditions
        c4 = np.concatenate(([riv], csp[:-1])) # river
        c1 = np.concatenate((cdp[1:], [ocn])) # ocean
        cc4 = np.concatenate(([riv], ccsp[:-1])) # river
        cc1 = np.concatenate((ccdp[1:], [ocn])) # ocean
        # effect of net fluxes
        cs = csp + bs*(q4*c4 + a21*q1*cdp - a34*q4*csp - q2*csp)
        cd = cdp + bd*(q1*c1 - a21*q1*cdp + a34*q4*csp - q3*cdp)
        ccs = ccsp + bs*(q4*cc4 + a21*q1*ccdp - a34*q4*ccsp - q2*ccsp) + dt*csp/86400
        ccd = ccdp + bd*(q1*cc1 - a21*q1*ccdp + a34*q4*ccsp - q3*ccdp) + dt*cdp/86400
        csp = cs.copy()
        cdp = cd.copy()
        ccsp = ccs.copy()
        ccdp = ccd.copy()
    
    return cs, cd, ccs, ccd

def plot_it(row=0, ylab='Tracer', xlab=False, leg=False, exp='blank'):
    ax = axes[row,0]
    ax.plot(xm/1e3, csa.T, '-c', linewidth = 1, alpha=.3)
    ax.plot(xm/1e3, cda.T, '-m', linewidth = 1, alpha=.3)
    ax.plot(xm/1e3, csa[-1,:], '-b', linewidth = 2, label = 'Surface')
    ax.plot(xm/1e3, cda[-1,:], '-r', linewidth = 2, label = 'Deep')
    ax.grid(True)
    ax.set_ylabel(ylab)
    ax.set_xlim(0,L/1e3)
    ax.set_ylim(bottom=0)
    if leg:
        ax.legend()
    if xlab:
        ax.set_xlabel('X (km)')
    ax.text(.95, .9, exp,
        transform=ax.transAxes, horizontalalignment='right')
    #
    ax = axes[row,1]
    ax.plot(td_vec, c_tot/V, '-c', linewidth=3)
    ax.plot(td_vec, c_tot_from_f/V, '-k')
    ax.set_xlim(0,ndays)
    ax.set_ylim(bottom=0)
    ax.grid(True)
    if xlab:
        ax.set_xlabel('Time (days)')
    ax.text(.95, .2, 'Volume Averaged '+ylab+' vs. Time',
        transform=ax.transAxes, horizontalalignment='right')
    # residence time
    c0 = c_tot[0]
    c1 = c_tot[-1]
    dc = np.abs(c1-c0)
    itres = np.argmax(np.abs(c_tot-c0) >= dc*(1-1/np.e))
    tres = td_vec[itres]
    try:
        ax.text(.95, .1, ('Residence Time = %0.1f days' % (tres)),
            transform=ax.transAxes, horizontalalignment='right')
    except OverflowError:
        pass
        
# Running the calculation function for different cases, and plotting
plt.close('all')

if False:
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(12,8))
    # a tracer that mimics salinity, as a test of the method
    # doubles as an experiment of how long it takes to get to the
    # steady state from zero
    csa, cda, td_vec, c_tot, c_tot_from_f = c_calc(ocn=Sin[-1])
    ax = axes[0,0]
    ax.plot(x/1e3, Sout, ':b', linewidth=3, label='Sout')
    ax.plot(x/1e3, Sin, '--r', linewidth=3, label='Sin')
    plot_it(row=0, ylab='Salinity', leg=True, exp='Reproduce Salinity', xlab=True)
    #
    ax = axes[1,0]
    wup = 86400*q1*a21/(B*dx)
    wdn = -86400*q4*a34/(B*dx)
    ax.plot(xm/1e3, wup, '-m', label='Wup (m/day)')
    ax.plot(xm/1e3, -wdn, '-g', label='-Wdown (m/day)')
    ax.plot(xm/1e3, wup + wdn, '-b', label='Wnet (m/day)')
    ax.plot(xm/1e3, wsink*np.ones(NX-1), '-c', label='-Wsink (m/day)')
    ax.legend()
    ax.grid(True)
    ax.set_xlim(0,L/1e3)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('X (km)')
    #
    ax = axes[1,1]
    ax.plot(x/1e3, Qin, '-r', linewidth=3, label='Qin (m3/s)')
    ax.plot(x/1e3, -Qout, '-b', linewidth=3, label='-Qout (m3/s)')
    ax.legend()
    ax.grid(True)
    ax.set_xlim(0,L/1e3)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('X (km)')
    # alternate residence time based on net c / (influx of c)
    tres_alt = c_tot[-1] / (Qin[-1]*socn*86400)
    ax.text(.95, .1, ('Alt. Residence Time = %0.1f days' % (tres_alt)),
        transform=ax.transAxes, horizontalalignment='right')
    ax.text(.95, .2, ('Freshwater Residence Time = %0.1f days' % (tres_fw)),
        transform=ax.transAxes, horizontalalignment='right')
    #
    # and an experiment on how long it takes to get rid of a salt-like tracer
    exp='salt_ic'
    csa, cda, td_vec, c_tot, c_tot_from_f = c_calc(ocn=0)
    plot_it(row=2, ylab='Salinity', leg=True, exp=exp, xlab=True)

# age
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12,8), squeeze=False)

cs, cd, ccs, ccd = age_calc(riv=1)
ax = axes[0,0]
Xm = xm/1e3
ax.plot(Xm, cs, label='cs')
ax.plot(Xm, cd, label='cd')
ax.legend()
ax.set_title('River Tracer')
ax.set_ylim(0,1)
ax.grid(True)
ax.set_xlim(0,L/1e3)

age_s = ccs/cs
age_d = ccd/cd
ax = axes[0,1]
ax.plot(Xm, age_s, label='age shallow')
ax.plot(Xm, age_d, label='age deep')
ax.legend()
ax.set_ylim(bottom=0)
ax.grid(True)
ax.set_xlim(0,L/1e3)

cs, cd, ccs, ccd = age_calc(ocn=1)
ax = axes[1,0]
Xm = xm/1e3
ax.plot(Xm, cs, label='cs')
ax.plot(Xm, cd, label='cd')
ax.legend()
ax.set_title('Ocean Tracer')
ax.set_ylim(0,1)
ax.grid(True)
ax.set_xlim(0,L/1e3)

age_s = ccs/cs
age_d = ccd/cd
ax = axes[1,1]
ax.plot(Xm, age_s, label='age shallow')
ax.plot(Xm, age_d, label='age deep')
ax.legend()
ax.set_ylim(bottom=0)
ax.grid(True)
ax.set_xlim(0,L/1e3)


if False:
    exp = ''
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(11,9))
    # a new tracer: river=1, ocean=0, with SINKING
    csa, cda, td_vec, c_tot, c_tot_from_f = c_calc(riv=1, sfact=sfact)
    plot_it(row=0, ylab='Tracer', exp='Tracer Source in River')
    # a new tracer: river=0, ocean=1, with SINKING
    csa, cda, td_vec, c_tot, c_tot_from_f = c_calc(ocn=1, sfact=sfact)
    plot_it(row=1, ylab='Tracer', xlab=True, exp='Tracer Source in Ocean')
    
plt.show()
