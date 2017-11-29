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
    params['Socn'] = 32
    params['Ssfc'] = 31
    params['k1'] = 1/(4*48)
    params['k2'] = 6.6/(80*48)
    return dims, params
    
    
def fSQ(Soutm, Sinp, Qr, Ut, L, H, A, params):
    pr = Bunch(params)
    c = np.sqrt(pr.g*pr.beta*pr.Socn*H)
    K = 0.028 * pr.Cd * Ut * H
    T = H**2 / K
    # find G = dsdx/Socn
    aa = pr.k2*c**2*T**2
    G = (-L + np.sqrt(L**2 + 8*aa*(Sinp-Soutm)/pr.Socn))/(4*aa)
    dq = A*pr.k1*c**2*T*G
    Qin = dq
    Qout = -(Qin + Qr)
    # adjust salinities to conserve salt flux
    ds = pr.Socn*aa*G**2
    eps = Qr*(Sinp-Soutm-2*ds)/(2*Qin+Qr)
    Sinm = Soutm + 2*ds - eps
    Soutp = Sinp - 2*ds - eps
    return Sinm, Soutp, Qin, Qout
    
def advance_s(Qr, Qr_p, Ut, Ut_p, dims, params, S1, S2, S1_p, S2_p, dt, do_check=False, landward_basin=True):
    
    # make easier access to dicts
    dm = Bunch(dims)
    pr = Bunch(params)
    
    # landward basin
    Soutm_p = S1_p
    Sinp_p = S2
    
    Sinm_p, Soutp_p, Qin_p, Qout_p = fSQ(Soutm_p, Sinp_p, Qr_p, Ut_p, dm.L_p, dm.H_p, dm.A_p, params)
    
    S1_pn = S1_p + dt*(S1_p*Qout_p + S2_p*Qin_p)/dm.V1_p
    S2_pn = S2_p + dt*(Sinm_p*Qin_p - S2_p*Qin_p)/dm.V2_p
    
    if landward_basin == False:
        S1_pn = 0
        S2_pn = 0
        Soutp_p = 0
        Qout_p = -Qr
        Qin_p = 0
        
    # seaward basin
    Soutm = S1
    Sinp = pr.Socn
    Qrr = Qr + Qr_p
    Sinm, Soutp, Qin, Qout = fSQ(Soutm, Sinp, Qrr, Ut, dm.L, dm.H, dm.A, params)
        
    Q21 = Qin-Qin_p
    if Q21 >= 0:
        F21 = Q21*S2
    else:
        F21 = Q21*S1
    S1n = S1 + dt*(S1*Qout - Soutp_p*Qout_p + F21)/dm.V1
    S2n = S2 + dt*(Sinm*Qin - S2*Qin_p  - F21)/dm.V2
    
    if do_check:
        # Checking on salt flux conservation at both ends
        # of the seaward strait
        # RESULT: it works perfectly now.
        Fl = Qin*Sinm + Qout*Soutm
        Fs = Qin*Sinp + Qout*Soutp
        print('F landward end = ' + str(Fl))
        print('F seaward end = ' + str(Fs))
        
    return S1n, S2n, S1_pn, S2_pn, Qin, Qin_p

def cubic_solver(a,b,c,d):
    """
    This gives the real root of a cubic polynomial of the form
    a*x**3 + b*x**2 + c*x + d = 0
    Parker MacCready 12/5/2002, recoded in python 10/25/2017
    """
    b3a = b/(3*a)
    b3a2 = b3a * b3a
    b3a3 = b3a * b3a * b3a
    p = 3*b3a2 - 2*(b/a)*b3a + c/a
    q = -b3a3 + (b/a)*b3a2 - (c/a)*b3a + d/a
    D = (p/3)**3 + (q/2)**2
    if D > 0:
        # calculate u and v (use the part of cube root)
        u_inn = -q/2 + np.sqrt(D)
        u = np.sign(u_inn)*np.abs(complex(u_inn)**(1/3))
        v_inn = -q/2 - np.sqrt(D)
        v = np.sign(v_inn)*np.abs(complex(v_inn)**(1/3))
    else:
        # calculate u and v (use an imaginary part of the cube root)
        u_inn = -q/2 + np.sqrt(complex(D))
        u = complex(u_inn)**(1/3)
        v_inn = -q/2 - np.sqrt(complex(D))
        v = complex(v_inn)**(1/3)
    # there are three roots, but we return only the first
    F1 = u + v 
    #F2 = -(u + v)/2 + (u - v)*i*np.sqrt(3)/2
    #F3 = -(u + v)/2 - (u - v)*i*np.sqrt(3)/2
    the_root = np.real(F1 - b3a)
    return the_root
    
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

