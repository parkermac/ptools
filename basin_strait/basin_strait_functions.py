import numpy as np
import matplotlib.pyplot as plt

class Bunch(object):
    # Class for streamlining dict access
    # Usage: x = Bunch(a_dict)
    # Then x.a is the same as a_dict['a']
    def __init__(self, a_dict):
        self.__dict__.update(a_dict)

def get_dims():
    # DIMENSIONS
    dims = dict()
    # seaward sill
    dims['H'] = 56 # depth
    dims['B'] = 9.7e3 # width
    dims['L'] = 35e3 # length
    dims['A'] = dims['B'] * dims['H']
    # seaward basin
    dims['V1'] = 77e9 * 0.2 # volume of upper basin box
    dims['V2'] = 77e9 - dims['V1'] # volume of deeper basin box
    # landward sill
    dims['H_p'] = 49 # depth
    dims['B_p'] = 1.4e3 # width
    dims['L_p'] = 9e3 # length
    dims['A_p'] = dims['B_p'] * dims['H_p']
    # landward basin
    dims['V1_p'] = 17e9 * 0.3 # volume of upper basin box
    dims['V2_p'] = 17e9 - dims['V1_p'] # volume of deeper basin box
    #
    # PARAMETERS
    params = dict()
    params['beta'] = 7.7e-4
    params['g'] = 9.8
    params['Cd'] = 2.6e-3
    # params['Socn'] = 32
    # params['Ssfc'] = 31
    params['k1'] = 1/(4*48)
    params['k2'] = 6.6/(80*48)
    return dims, params
    
    
def fSQ(Stopm, Sbotp, Qr, Ut, L, H, A, params, Socn):
    pr = Bunch(params)
    c = np.sqrt(pr.g*pr.beta*Socn*H)
    K = 0.028 * pr.Cd * Ut * H
    T = H**2 / K
    # find G = dsdx/Socn
    aa = pr.k2*c**2*T**2
    G = (-L + np.sqrt(L**2 + 8*aa*(Sbotp-Stopm)/Socn))/(4*aa)
    dq = A*pr.k1*c**2*T*G
    Qbot = dq
    Qtop = -dq - Qr
    # adjust salinities to conserve salt flux
    ds = Socn*aa*G**2
    eps = Qr*(Sbotp-Stopm-2*ds)/(2*Qbot+Qr)
    Sbotm = Stopm + 2*ds - eps
    Stopp = Sbotp - 2*ds - eps
    return Sbotm, Stopp, Qbot, Qtop
    
def fSQ_rev(Stopp, Sbotm, Qr, Ut, L, H, A, params, Socn):
    # for the case of reversed flow
    pr = Bunch(params)
    c = np.sqrt(pr.g*pr.beta*Socn*H)
    K = 0.028 * pr.Cd * Ut * H
    T = H**2 / K
    # find G = dsdx/Socn
    aa = pr.k2*c**2*T**2
    G = -(-L + np.sqrt(L**2 + 8*aa*(Sbotm-Stopp)/Socn))/(4*aa)
    # note that G is now negative, so dq is negative as well
    dq = A*pr.k1*c**2*T*G
    Qbot = dq # so this is negative: a sink from the basin
    Qtop = -dq - Qr
    # adjust salinities to conserve salt flux
    ds = Socn*aa*G**2
    eps = Qr*(Sbotm-Stopp-2*ds)/(2*Qbot+Qr)
    Sbotp = Stopp + 2*ds - eps
    Stopm = Sbotm - 2*ds - eps
    return Sbotp, Stopm, Qbot, Qtop
    
def advance_s(Qr, Qr_p, Ut, Ut_p, Socn, Ssfc,
        dims, params, S1, S2, S1_p, S2_p, dt,
        do_check=False, landward_basin=True):
    
    # make easier access to dicts
    dm = Bunch(dims)
    pr = Bunch(params)
    
    # landward basin
    Stopm_p = S1_p
    Sbotp_p = S2
    Stopp_p = S1
    Sbotm_p = S2_p
    # which gradient is biggest?
    gr_norm_p = Sbotp_p - Stopm_p
    gr_rev_p = Sbotm_p - Stopp_p
    rev_p = False
    ff = 2
    # if gr_norm_p >= ff*gr_rev_p:
    #     rev_p = False
    if gr_rev_p > ff*gr_norm_p:
        rev_p = True
        
    if landward_basin == False:
        S1_pn = 0
        S2_pn = 0
        Stopp_p = 0
        Qtop_p = -Qr
        Qbot_p = 0
    else:
        if rev_p == False:
            # Normal case
            Sbotm_p, Stopp_p, Qbot_p, Qtop_p = fSQ(Stopm_p, Sbotp_p, Qr_p, Ut_p,
                dm.L_p, dm.H_p, dm.A_p, params, Socn)
            S1_pn = S1_p + dt*(S1_p*Qtop_p + S2_p*Qbot_p)/dm.V1_p
            S2_pn = S2_p + dt*(Sbotm_p*Qbot_p - S2_p*Qbot_p)/dm.V2_p
        elif rev_p == True:
            # Reversed case
            Stopp_p = S1
            Sbotm_p = S2_p
            Sbotp_p, Stopm_p, Qbot_p, Qtop_p = fSQ_rev(Stopp_p, Sbotm_p, Qr_p, Ut_p,
                dm.L_p, dm.H_p, dm.A_p, params, Socn)
            S1_pn = S1_p + dt*(Stopm_p*Qtop_p + S1_p*Qbot_p)/dm.V1_p
            S2_pn = S2_p + dt*(S2_p*Qbot_p - S1_p*Qbot_p)/dm.V2_p
        
    # seaward basin
    Stopm = S1
    Sbotp = Socn
    Qrr = Qr + Qr_p
    
    # landward basin
    Stopm = S1
    Sbotp = Socn
    Stopp = Ssfc
    Sbotm = S2
    # which gradient is biggest?
    gr_norm = Sbotp - Stopm
    gr_rev = Sbotm - Stopp
    if gr_norm >= gr_rev:
        rev = False
    elif gr_rev > gr_norm:
        rev = True
    
    if rev == False:
        # Normal case
        Sbotm, Stopp, Qbot, Qtop = fSQ(Stopm, Sbotp, Qrr, Ut, dm.L, dm.H, dm.A, params, Socn)
        Q21 = Qbot-Qbot_p
        if Q21 >= 0:
            F21 = Q21*S2
        else:
            F21 = Q21*S1
        S1n = S1 + dt*(S1*Qtop - Stopp_p*Qtop_p + F21)/dm.V1
        S2n = S2 + dt*(Sbotm*Qbot - S2*Qbot_p  - F21)/dm.V2
    elif rev == True:
        # Reversed case
        Stopp = Ssfc
        Sbotm = S2
        Sbotp, Stopm, Qbot, Qtop = fSQ_rev(Stopp, Sbotm, Qrr, Ut, dm.L, dm.H, dm.A, params, Socn)
        Q21 = Qbot-Qbot_p
        if Q21 >= 0:
            F21 = Q21*S2
        else:
            F21 = Q21*S1
        if rev_p == True:
            S1n = S1 + dt*(Stopm*Qtop - S1*Qtop_p + F21)/dm.V1
            S2n = S2 + dt*(S2*Qbot - Sbotp_p*Qbot_p  - F21)/dm.V2
        elif rev_p == False:
            S1n = S1 + dt*(Stopm*Qtop - Stopp_p*Qtop_p + F21)/dm.V1
            S2n = S2 + dt*(S2*Qbot - S2*Qbot_p  - F21)/dm.V2
    
    if do_check:
        # Checking on salt flux conservation at both ends
        # of the seaward strait
        # RESULT: it works perfectly now.
        Fl = Qbot*Sbotm + Qtop*Stopm
        Fs = Qbot*Sbotp + Qtop*Stopp
        print('F landward end = ' + str(Fl))
        print('F seaward end = ' + str(Fs))
        
    return S1n, S2n, S1_pn, S2_pn, Qbot, Qbot_p
    
# RC SETUP (plotting defaults)
def set_rc(fs_big, fs_small, lw_big, lw_small):
    plt.rc('xtick', labelsize=fs_small)
    plt.rc('ytick', labelsize=fs_small)
    plt.rc('xtick.major', size=10, pad=5, width=lw_small)
    plt.rc('ytick.major', size=10, pad=5, width=lw_small)
    plt.rc('axes', lw=lw_small)
    plt.rc('lines', lw=lw_big)
    plt.rc('font', size=fs_big)
    plt.rc('grid', color='g', ls='-', lw=lw_small, alpha=.3)
    
def add_line(ax, nd):
    aa = ax.axis()
    ax.plot([nd, nd], aa[-2:], '-k')

if __name__ == '__main__':
    # example use of the functions
    x = cubic_solver(1,2,3,4)
    print('The real root of the cubic is ' + str(the_root))

