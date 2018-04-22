"""
Functions for evolution of a system of basins and straits.

Both the basins and the straits have two layers, 1 = upper
and 2 = lower.
"""

import numpy as np
from collections import OrderedDict as od
import pandas as pd

def make_BaSt(case, names):
    """
    This function is where we set the basin and strait dimensions
    and their arrangement of connections.
    
    Ba and St are ordered dicts of ordered dicts.  The  primary keys
    are basin keys such as 'b' or strait keys such as 'ba'.
    The strait keys identify which two basins the strait connects,
    ordered landward-seaward.
    
    The actual creation of the dicts-of-dicts is done by the
    functions make_Ba() and make_St().
    
    What we pass to these functions are the temporary dicts bbdict
    and ssdict, which are just compact ways of initially organizing
    the basin and strait dimensions to make it easy to create cases.
    
    The assumed packing order of the dictionary entries are
    bbdict = {basin key: [V1 (m3), V2 (m3)], ...}
    ssdict = {strait key: [H (m), B (m), L (m), r (string)], ...}
    
    The 'r' string in each ssdict entry is all the basins landward of
    a strait, from which we compute a net river flow.
    """
    bbdict = od()
    ssdict = od()
    if case == 'base':
        bbdict = {'a':[1e12, 2e12],
                'b':[10e9, 70e9],
                'c':[2e9, 14e9],
                'd':[2e9, 14e9]}
        ssdict = {'ba': [30, 3e3, 10e3, 'bcd'],
                'cb': [30, 3e3, 10e3, 'cd'],
                'dc': [30, 3e3, 10e3, 'd']}
    # make the initial dicts-of-dicts
    Ba = make_Ba(bbdict)
    St = make_St(ssdict)
    # then add variables
    for name in names:
        if name == 's':
            v = pr['Socn']
        else:
            v = 0
        NAME = name.upper()
        for bk in Ba.keys():
            binfo = Ba[bk]
            binfo[NAME+'1'] = v
            binfo[NAME+'2'] = v
            Ba[bk] = binfo
        name = name.lower()
        for sk in St.keys():
            sinfo = St[sk]
            tlist = ['1m', '1p', '2m', '2p']
            for t in tlist:
                sinfo[name+t] = v
            St[sk] = sinfo
    return Ba, St

def set_Qr(Ba, t, fcase='none'):
    """
    This sets the river flow entering each basin, as
    a function of time.
    """
    td = t/86400
    if fcase == 'base':
        for bk in Ba.keys():
            if bk=='c':
                Ba[bk]['Qr'] = 500 + 2000*np.exp(-((td-500)/30)**2)
            else:
                Ba[bk]['Qr'] = 500 + 500*np.exp(-((td-500)/30)**2)
    return Ba

def set_Ut(St, t, fcase='none'):
    """
    This sets the tidal velocity in each strait, as
    a function of time.
    """
    td = t/86400
    if fcase == 'base':
        for sk in St.keys():
            St[sk]['Ut'] = 1 + 0.1*np.sin(2*np.pi*td/14)
    return St

def get_pr():
    """
    This sets parameters that are shared by various calculations,
    especially in the straits.
    """
    pr = dict()
    pr['g'] = 9.8; pr['beta'] = 7.7e-4; pr['Cd'] = 2.6e-3
    pr['k1'] = 1/(4*48); pr['k2'] = 6.6/(80*48);
    pr['Socn'] = 30; pr['Sriv'] = 0
    pr['dt'] = int(86400) # integration time step (sec)
    return pr
pr = get_pr() # make pr available to the other functions here

def make_Ba(bbdict):
    """
    This is where we create the dict-of-dicts that the code uses,
    for the basins.  We also add some other state variables with
    default initial values.
    """
    Ba = od()
    for bk in bbdict.keys():
        binfo = dict()
        binfo['V1'] = bbdict[bk][0]
        binfo['V2'] = bbdict[bk][1]
        binfo['Qr'] = 500
        Ba[bk] = binfo
    return Ba

def make_St(ssdict):
    """
    This is where we create the dict-of-dicts that the code uses,
    for the straits.  We also add some other state variables with
    default initial values.
    """
    St = od()
    for sk in ssdict.keys():
        sinfo = dict()
        sinfo['H'] = ssdict[sk][0]
        sinfo['B'] = ssdict[sk][1]
        sinfo['L'] = ssdict[sk][2]
        sinfo['r'] = ssdict[sk][3]
        sinfo['Ut'] = 1
        # the transports have one value in each layer
        sinfo['q1'] = 0
        sinfo['q2'] = 0
        sinfo['rev'] = False
        St[sk] = sinfo
    return St
    
def make_BaSt_out(Ba, St, NT, names):
    """
    Initialize the output arrays, and make dicts to specify
    which variable goes in which column.
    """
    v = ['Qr']
    for name in names:
        NAME = name.upper()
        v.append(NAME+'1')
        v.append(NAME+'2')
    nv = len(v)
    cc = []
    vv = []
    for k in Ba.keys():
        for ii in range(nv):
            cc.append(k)
        vv += v
    Ba_out = np.nan * np.ones((NT, len(cc)))
    cv = []
    for ii in range(len(cc)):
        cv.append((cc[ii],vv[ii]))
    Ba_out_cols = dict(zip(cv, range(len(cc))))

    v = ['Ut', 'q1', 'q2']
    for name in names:
        name = name.lower()
        tlist = ['1m', '1p', '2m', '2p']
        for t in tlist:
            v.append(name+t)
    nv = len(v)
    cc = []
    vv = []
    for k in St.keys():
        for ii in range(nv):
            cc.append(k)
        vv += v
    St_out = np.nan * np.ones((NT, len(cc)))
    cv = []
    for ii in range(len(cc)):
        cv.append((cc[ii],vv[ii]))
    St_out_cols = dict(zip(cv, range(len(cc))))

    return Ba_out, Ba_out_cols, St_out, St_out_cols
    
def strait(Qr_sum, sinfo, pr):
    """
    This updates the four salinity values and two transports
    that define the dyanmical state of each of the straits.
    """
    g = pr['g']
    beta = pr['beta']
    Socn = pr['Socn']
    Cd = pr['Cd']
    k1 = pr['k1']
    k2 = pr['k2']
    H = sinfo['H']
    A = sinfo['H'] * sinfo['B']
    L = sinfo['L']
    Ut = sinfo['Ut']
    prev_rev = sinfo['rev']
    c = np.sqrt(g * beta * Socn * H)
    K = 0.028 * Cd * Ut * H
    T = H**2 / K
    # which gradient is biggest?
    gr_norm = sinfo['s2p'] - sinfo['s1m']
    gr_rev = sinfo['s2m'] - sinfo['s1p']
    ff = 1.5 # factor to damp oscillations in the straits
    rev = prev_rev
    # if gr_norm_p >= ff*gr_rev_p:
    if (gr_rev > ff*gr_norm) and prev_rev==False :
        rev = True
    elif (gr_norm > ff*gr_rev) and prev_rev==True :
        rev = False
    # find G = dsdx/Socn
    if rev==False:
        aa = k2 * c**2 * T**2
        G = (-L + np.sqrt(L**2 + 8*aa*(sinfo['s2p']-sinfo['s1m'])/Socn))/(4*aa)
        dq = A*k1*c**2*T*G
        # sign convention for fluxes: positive seaward
        q2 = -dq
        q1 = dq + Qr_sum
        # adjust salinities to conserve salt flux
        ds = Socn*aa*G**2
        eps = Qr_sum*(sinfo['s2p']-sinfo['s1m']-2*ds)/(-2*q2+Qr_sum)
        sinfo['s2m'] = sinfo['s1m'] + 2*ds - eps
        sinfo['s1p'] = sinfo['s2p'] - 2*ds - eps
    elif rev==True:
        aa = k2 * c**2 * T**2
        G = -(-L + np.sqrt(L**2 + 8*aa*(sinfo['s2m']-sinfo['s1p'])/Socn))/(4*aa)
        dq = A*k1*c**2*T*G
        q2 = -dq
        q1 = dq + Qr_sum
        # adjust salinities to conserve salt flux
        ds = Socn*aa*G**2
        eps = Qr_sum*(sinfo['s2m']-sinfo['s1p']-2*ds)/(-2*q2+Qr_sum)
        sinfo['s2p'] = sinfo['s1p'] + 2*ds - eps
        sinfo['s1m'] = sinfo['s2m'] - 2*ds - eps
    sinfo['q1'] = q1
    sinfo['q2'] = q2
    sinfo['rev'] = rev
    return sinfo

def bcalc(Ba, St, Ba_out, Ba_out_cols, St_out, St_out_cols, NT, fcase, ic=False):
    """
    This is the main calculate that integrates the state
    of all the basins forward over NT timesteps.
    
    """
    sv_list = [k for (x,k) in St_out_cols.keys() if x=='ba']
    bv_list = [k for (x,k) in Ba_out_cols.keys() if x=='a']
    pr = get_pr()
    TD = np.nan * np.ones(NT)
    # time integration
    dt = pr['dt']
    t = 0
    ii = 0
    while ii < NT:
        # update forcing
        if ic == True:
            Ba = set_Qr(Ba, 0, fcase=fcase)
            St = set_Ut(St, 0, fcase=fcase)
        else:
            Ba = set_Qr(Ba, t, fcase=fcase)
            St = set_Ut(St, t, fcase=fcase)
        # first update all the straits
        for sk in St.keys():
            sinfo = St[sk]
            Qr_sum = 0
            for rbk in sinfo['r']:
                Qr_sum += Ba[rbk]['Qr']
            sinfo = strait(Qr_sum, sinfo, pr)
            St[sk] = sinfo
            # save results
            for v in sv_list:
                St_out[ii, St_out_cols[(sk,v)]] = sinfo[v]
        # go through all the basins
        for bk in Ba.keys():
            binfo = Ba[bk]
            V1 = binfo['V1']
            V2 = binfo['V2']
            # save results
            for v in bv_list:
                Ba_out[ii, Ba_out_cols[(bk,v)]] = binfo[v]
            TD[ii] = t/86400
            if bk == 'a':
                pass
            else:
                # find all fluxes from straits
                qin1 = 0
                sqin1 = binfo['Qr'] * pr['Sriv']
                # add river salt flux for generality
                qin2 = 0
                sqin2 = 0
                for sk in St.keys():
                    sinfo = St[sk]
                    rev = sinfo['rev']
                    if sk[0] == bk: # strait is seaward of basin
                        if rev == False:
                            qin1 += -sinfo['q1']
                            sqin1 += -sinfo['q1']*binfo['S1']
                            qin2 += -sinfo['q2']
                            sqin2 += -sinfo['q2']*sinfo['s2m']
                        elif rev == True:
                            qin1 += -sinfo['q1']
                            sqin1 += -sinfo['q1']*sinfo['s1m']
                            qin2 += -sinfo['q2']
                            sqin2 += -sinfo['q2']*binfo['S2']
                    elif sk[1] == bk: # strait is landward of basin
                        if rev == False:
                            qin1 += sinfo['q1']
                            sqin1 += sinfo['q1']*sinfo['s1p']
                            qin2 += sinfo['q2']
                            sqin2 += sinfo['q2']*binfo['S2']
                        elif rev == True:
                            qin1 += sinfo['q1']
                            sqin1 += sinfo['q1']*binfo['S1']
                            qin2 += sinfo['q2']
                            sqin2 += sinfo['q2']*sinfo['s2p']
                    else: # strait is not connected to basin
                        pass
                # calculate vertical salt transport in basin
                #print('%s: qin1 = %d, qin2 = %d' % (bk, qin1, qin2))
                if qin2 > 0:
                    squp = qin2 * binfo['S2']
                elif qin2 <= 0:
                    squp = qin2 * binfo['S1']
                # update basin salinity
                Ba[bk]['S1'] = binfo['S1'] + dt*(sqin1 + squp)/V1
                Ba[bk]['S2'] = binfo['S2'] + dt*(sqin2 - squp)/V2
        # reset upstream values for all straits
        for sk in St.keys():
            sinfo = St[sk]
            sinfo['s1m'] = Ba[sk[0]]['S1']
            sinfo['s2m'] = Ba[sk[0]]['S2']
            sinfo['s1p'] = Ba[sk[1]]['S1']
            sinfo['s2p'] = Ba[sk[1]]['S2']
            St[sk] = sinfo
        t += dt
        ii += 1
    return (Ba, St, Ba_out, St_out, TD)

