"""
Functions for the rflx code.
"""

import numpy as np

# function to create Sin and Sout
def get_Sio_chatwin(Qr, Socn, ds, nx):
    L0 = 50e3 # length of salt intrusion for Qr = Qr0
    Qr0 = 1e3 # m3/s
    L= L0 * (Qr/Qr0)**(-1/3) # length of salt intrusion (m)
    a = Socn/(L**1.5)
    alpha = ds/L
    x = np.linspace((alpha/(2*a))**2,L,nx)
    Sin = a*x**1.5 + alpha*x/2
    Sout = a*x**1.5 - alpha*x/2
    return Sin, Sout, x, L

def a_calc(Sin, Sout):
    Sin1 = Sin[1:]
    Sin0 = Sin[:-1]
    Sout1 = Sout[1:]
    Sout0 = Sout[:-1]
    # Calculate Efflux-Reflux fractions (a1 and a0), at xm
    # (notation different from Cokelet and Stewart 1985)
    a0 = (Sout0/Sin0)*(Sin1-Sin0)/(Sin1-Sout0) # reflux - down
    a1 = (Sin1/Sout1)*(Sout1-Sout0)/(Sin1-Sout0) # efflux - up
    return a0, a1
    
def get_time_step(dx, Qout, B, hs, ndays):
    """
    Calculate the time step dynamically using some factor (0.9) times
    the shortest time for water to traverse a grid box in the
    along_channel direction (experiments showed that it was unstable
    for cfl_factor > 1.1).
    """
    cfl_factor = 0.9
    dt = cfl_factor * np.min(dx/(Qout[1:]/(B*hs))) # time step (s)
    NT = int(ndays*86400/dt) # total number of time steps
    NS = 10 # number of saves
    return dt, NT, NS

def c_calc(csp, cdp, info_tup, riv=0, ocn=0, Ts=np.inf, Td=np.inf, Cs=1, Cd=0):
    """
    This is the main computational engine for the time dependent solution.
    
    csp, cdp = vectors [xm] of initial tracer concentrations in the two layers
    csa, cda = arrays [time, xm] of the surface and deep concentrations
        over the course of the simulation
    
    info_tup = tuple of properties defining the grid, timestep and circulation
    
    riv = value of tracer coming in from river
    ocn = value of tracer coming in from ocean
    Ts = relaxation timescale [days] for surface layer
    Td = relaxation timescale [days] for deep layer
    Cs = value of tracer to relax to in surface layer
    Cd = value of tracer to relax to in deep layer
    """
    # unpack some parameters
    NS, NX, NT, dt, dvs, dvd, Qout, Qin, a0, a1 = info_tup
    # Get Q's at box edges used in the box model.
    # NOTE: all Q's are positive for the box model, whereas
    # Qout is negative in my usual TEF sign convention.
    Qout0 = Qout[:-1]
    Qin0 = Qin[:-1]
    Qout1 = Qout[1:]
    Qin1 = Qin[1:]
    NTs = int(NT/NS) # save interval, in timesteps
    # initialize arrays to save in
    csa = np.nan * np.ones((NS,NX-1))
    cda = np.nan * np.ones((NS,NX-1))
    tta = 0 # index for periodic saves
    for tt in range(NT):
        # boundary conditions
        Cout0 = np.concatenate(([riv], csp[:-1])) # river
        Cin1 = np.concatenate((cdp[1:], [ocn])) # ocean
        # update fields
        cs = csp + (dt/dvs)*( Qout0*(1-a0)*Cout0 + Qin1*a1*Cin1 - Qout1*csp) + (Cs-csp)*(dt/86400)/Ts
        cd = cdp + (dt/dvd)*( Qin1*(1-a1)*Cin1 + Qout0*a0*Cout0 - Qin0*cdp) + (Cd-cdp)*(dt/86400)/Td
        cs[cs<0] = 0
        cd[cd<0] = 0
        cd[0] = cd[1] # this helps when using riv = ocn = const.
        csp = cs.copy()
        cdp = cd.copy()
        if (np.mod(tt, NTs) == 0) and tta < NS:
            # periodic save
            csa[tta,:] = cs
            cda[tta,:] = cd
            tta += 1
    return csa, cda
    
