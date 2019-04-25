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

def qa_calc(Qin, Qout, Sin, Sout, I):
    # Use I = 0 to ignore the second river part
    
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

    # add the upwelling required for the second river
    if I > 0:
        a21x = np.zeros_like(a21)
        a21x[I] = (q1[I]-q1[I-1])/q1[I]
        a21 = a21 + a21x
    
    return q1, q2, q3, q4, a21, a34

def c_calc(csp, cdp, info_tup, riv=0, ocn=0, riv2=0, do_age=False):
    # unpack some parameters
    NS, NX, NT, dt, dvs, dvd, q1, q2, q3, q4, a21, a34, o_Qr, I = info_tup
    NTs = int(NT/NS) # save interval, in timesteps
    if do_age:
        # arrays for age calculation
        ccsp = csp.copy()
        ccdp = cdp.copy()
    # arrays to save in
    csa = np.nan * np.ones((NS,NX-1))
    cda = np.nan * np.ones((NS,NX-1))
    # time series to save in
    c_tot = np.nan * np.ones(NT) # net amount of tracer
    t_vec = np.nan * np.ones(NT) # time (s)
    # these will be time series of the tracer flux at open boundaries
    f1_vec = np.nan * np.ones(NT)
    f2_vec = np.nan * np.ones(NT)
    f3_vec = np.nan * np.ones(NT)
    f4_vec = np.nan * np.ones(NT)
    friv2_vec = np.nan * np.ones(NT)
    #
    riv_mask = np.zeros(NX-1)
    riv_mask[I] = 1
        
    tta = 0 # index for periodic saves
    for tt in range(NT):
        # boundary conditions
        c4 = np.concatenate(([riv], csp[:-1])) # river
        c1 = np.concatenate((cdp[1:], [ocn])) # ocean
        if do_age:
            cc4 = np.concatenate(([riv], ccsp[:-1])) # river
            cc1 = np.concatenate((ccdp[1:], [ocn])) # ocean
        # save time series entries
        t_vec[tt] = tt*dt
        c_tot[tt] = np.sum(csp*dvs + cdp*dvd)
        f1_vec[tt] = q1[-1]*c1[-1]
        f2_vec[tt] = - q2[-1]*csp[-1]
        f3_vec[tt] = - q3[0]*cdp[0]
        f4_vec[tt] = q4[0]*c4[0]
        friv2_vec[tt] = o_Qr*riv2
        # update fields
        cs = csp + (dt/dvs)*(q4*c4 + a21*q1*cdp - a34*q4*csp - q2*csp + o_Qr*riv2*riv_mask)
        cd = cdp + (dt/dvs)*(q1*c1 - a21*q1*cdp + a34*q4*csp - q3*cdp)
        if do_age:
            # ageing versions
            ccs = ccsp + (dt/dvs)*(q4*cc4 + a21*q1*ccdp - a34*q4*ccsp - q2*ccsp + o_Qr*riv2*riv_mask) + dt*csp/86400
            ccd = ccdp + (dt/dvs)*(q1*cc1 - a21*q1*ccdp + a34*q4*ccsp - q3*ccdp) + dt*cdp/86400
            ccsp = ccs.copy()
            ccdp = ccd.copy()
        csp = cs.copy()
        cdp = cd.copy()
        if (np.mod(tt, NTs) == 0) and tta < NS:
            # periodic save
            csa[tta,:] = cs
            cda[tta,:] = cd
            tta += 1
    # work on time series
    T = t_vec/86400 # time axis in days
    # pack things
    f_tup = (T, c_tot, f1_vec, f2_vec, f3_vec, f4_vec, friv2_vec)
    
    if do_age:
        age_tup = (cs, cd, ccs, ccd)
        return csa, cda, f_tup, age_tup
    else:
        return csa, cda, f_tup