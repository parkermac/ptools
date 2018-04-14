"""
Functions for a system of basins and straits.
"""

import numpy as np
from collections import OrderedDict as od

def get_BaSt(case):
    # packing order
    # {basin key: [V1 (m3), V2 (m3)]}
    # {strait key: [H (m), B (m), L (m), r (string)]}
    bbdict = od()
    ssdict = od()
    if case == 'base':
        bbdict = {'a':[1e12, 2e12], 'b':[10e9, 70e9], 'c':[2e9, 14e9], 'd':[2e9, 14e9]}
        ssdict = {'ba': [30, 3e3, 10e3, 'bcd'], 'cb': [30, 3e3, 10e3, 'cd'], 'dc': [30, 3e3, 10e3, 'd']}
    Ba = make_Ba(bbdict)
    St = make_St(ssdict)
    return Ba, St
    
def set_Qr(Ba, t, fcase='none'):
    td = t/86400
    if fcase == 'base':
        for bk in Ba.keys():
            if bk=='c':
                Ba[bk]['Qr'] = 500 + 1000*np.exp(-((td-500)/60)**2)
            else:
                Ba[bk]['Qr'] = 500
    return Ba
    
def set_Ut(St, t, fcase='none'):
    td = t/86400
    if fcase == 'base':
        for sk in St.keys():
            St[sk]['Ut'] = 1 + 0.1*np.sin(2*np.pi*td/14)
    return St

def get_pr():
    # parameters
    pr = dict()
    pr['g'] = 9.8; pr['beta'] = 7.7e-4; pr['Cd'] = 2.6e-3
    pr['k1'] = 1/(4*48); pr['k2'] = 6.6/(80*48); pr['Socn'] = 30
    return pr
pr = get_pr()

def make_Ba(bbdict):
    Ba = od()
    for bk in bbdict.keys():
        bdict = dict()
        bdict['V1'] = bbdict[bk][0]
        bdict['V2'] = bbdict[bk][1]
        bdict['S1'] = pr['Socn']
        bdict['S2'] = pr['Socn']
        bdict['Qr'] = 500
        Ba[bk] = bdict
    return Ba
        
def make_St(ssdict):
    St = od()
    for sk in ssdict.keys():
        sdict = dict()
        sdict['H'] = ssdict[sk][0]
        sdict['B'] = ssdict[sk][1]
        sdict['L'] = ssdict[sk][2]
        sdict['r'] = ssdict[sk][3]
        sdict['Ut'] = 1
        sdict['s1'] = [pr['Socn'], pr['Socn']]
        sdict['s2'] = [pr['Socn'], pr['Socn']]
        sdict['q1'] = 0
        sdict['q2'] = 0
        sdict['rev'] = False
        St[sk] = sdict
    return St

def strait(Qr_sum, sinfo, pr):
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
    s1 = sinfo['s1']
    s2 = sinfo['s2']
    prev_rev = sinfo['rev']
    c = np.sqrt(g * beta * Socn * H)
    K = 0.028 * Cd * Ut * H
    T = H**2 / K
    # p "plus" is seaward
    # m "minus" is landward
    s1m = s1[0]
    s1p = s1[1]
    s2m = s2[0]
    s2p = s2[1]
    # which gradient is biggest?
    gr_norm = s2p - s1m
    gr_rev = s2m - s1p
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
        G = (-L + np.sqrt(L**2 + 8*aa*(s2p-s1m)/Socn))/(4*aa)
        dq = A*k1*c**2*T*G
        # sign convention for fluxes: positive seaward
        q2 = -dq
        q1 = dq + Qr_sum
        # adjust salinities to conserve salt flux
        ds = Socn*aa*G**2
        eps = Qr_sum*(s2p-s1m-2*ds)/(-2*q2+Qr_sum)
        s2m = s1m + 2*ds - eps
        s1p = s2p - 2*ds - eps
    elif rev==True:
        aa = k2 * c**2 * T**2
        G = -(-L + np.sqrt(L**2 + 8*aa*(s2m-s1p)/Socn))/(4*aa)
        dq = A*k1*c**2*T*G
        q2 = -dq
        q1 = dq + Qr_sum
        # adjust salinities to conserve salt flux
        ds = Socn*aa*G**2
        eps = Qr_sum*(s2m-s1p-2*ds)/(-2*q2+Qr_sum)
        s2p = s1p + 2*ds - eps
        s1m = s2m - 2*ds - eps
    sinfo['s1'] = [s1m, s1p]
    sinfo['s2'] = [s2m, s2p]
    sinfo['q1'] = q1
    sinfo['q2'] = q2
    sinfo['rev'] = rev
    return sinfo

def bcalc(St, Ba, NT, fcase, ic=False):
    pr = get_pr()
    # initialize output
    S11 = dict()
    S22 = dict()
    Qrr = dict()
    for bk in Ba.keys():
        S11[bk] = np.nan * np.ones(NT)
        S22[bk] = np.nan * np.ones(NT)
        Qrr[bk] = np.nan * np.ones(NT)
    Q11 = dict()
    Q22 = dict()
    Utt = dict()
    for sk in St.keys():
        Q11[sk] = np.nan * np.ones(NT)
        Q22[sk] = np.nan * np.ones(NT)
        Utt[sk] = np.nan * np.ones(NT)
    TD = np.nan * np.ones(NT)
    # time integration
    dt = 86400
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
            Q11[sk][ii] = sinfo['q1']
            Q22[sk][ii] = sinfo['q2']
            Utt[sk][ii] = sinfo['Ut']
        # go through all the basins
        for bk in Ba.keys():
            binfo = Ba[bk]
            V1 = binfo['V1']
            V2 = binfo['V2']
            S1 = binfo['S1']
            S2 = binfo['S2']
            # save results
            S11[bk][ii] = S1
            S22[bk][ii] = S2
            Qrr[bk][ii] = binfo['Qr']
            TD[ii] = t/86400
            if bk == 'a':
                pass
            else:
                # find all fluxes from straits
                qin1 = 0
                sqin1 = 0
                qin2 = 0
                sqin2 = 0
                for sk in St.keys():
                    sinfo = St[sk]
                    rev = sinfo['rev']
                    if sk[0] == bk: # strait is seaward of basin
                        if rev == False:
                            qin1 += -sinfo['q1']
                            sqin1 += -sinfo['q1']*S1
                            qin2 += -sinfo['q2']
                            sqin2 += -sinfo['q2']*sinfo['s2'][0]
                        elif rev == True:
                            qin1 += -sinfo['q1']
                            sqin1 += -sinfo['q1']*sinfo['s1'][0]
                            qin2 += -sinfo['q2']
                            sqin2 += -sinfo['q2']*S2
                    elif sk[1] == bk: # strait is landward of basin
                        if rev == False:
                            qin1 += sinfo['q1']
                            sqin1 += sinfo['q1']*sinfo['s1'][1]
                            qin2 += sinfo['q2']
                            sqin2 += sinfo['q2']*S2
                        elif rev == True:
                            qin1 += sinfo['q1']
                            sqin1 += sinfo['q1']*S1
                            qin2 += sinfo['q2']
                            sqin2 += sinfo['q2']*sinfo['s2'][1]
                    else: # strait is not connected to basin
                        pass
                # calculate vertical salt transport in basin
                #print('%s: qin1 = %d, qin2 = %d' % (bk, qin1, qin2))
                if qin2 > 0:
                    squp = qin2 * S2
                elif qin2 <= 0:
                    squp = qin2 * S1
                # update basin salinity
                Ba[bk]['S1'] = S1 + dt*(sqin1 + squp)/V1
                Ba[bk]['S2'] = S2 + dt*(sqin2 - squp)/V2
        # reset upstream values for all straits
        for sk in St.keys():
            sinfo = St[sk]
            sinfo['s1'] = [Ba[sk[0]]['S1'], Ba[sk[1]]['S1']]
            sinfo['s2'] = [Ba[sk[0]]['S2'], Ba[sk[1]]['S2']]
            St[sk] = sinfo
        t += dt
        ii += 1
    return (S11, S22, Qrr, Q11, Q22, Utt, TD, St, Ba)
    
def print_strait(sinfo):
    s1 = sinfo['s1']
    s2 = sinfo['s2']
    q1 = sinfo['q1']
    q2 = sinfo['q2']
    rev = sinfo['rev']
    print('Reversed = %s' % (str(rev)))
    print('s1m = %0.2f, s1p = %0.2f' % (s1[0], s1[1]))
    print('s2m = %0.2f, s2p = %0.2f' % (s2[0], s2[1]))
    print('q1 = %d, q2 = %d' % (int(q1), int(q2)))
    Fm = q1*s1[0] + q2*s2[0]
    Fp = q1*s1[1] + q2*s2[1]
    print('Fm = %0.1f, Fp = %0.1f' % (Fm, Fp))

def test_strait(Qr_sum, sinfo, pr):
    # testing strait()
    sinfo['s1'] = [10, 12]
    sinfo['s2'] = [15, 17]
    sinfo['q1'] = 0
    sinfo['q2'] = 0
    sinfo['rev'] = False
    print('\nBefore:')
    print_strait(sinfo)
    sinfo = strait(Qr_sum, sinfo, pr)
    print('\nAfter:')
    print_strait(sinfo)
