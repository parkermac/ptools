import numpy as np

def get_dims():
    # SIZES
    # seaward sill
    H = 56 # depth
    B = 9.7e3 # width
    L = 35e3 # length
    # seaward basin
    V1 = 77e9 * 0.2 # volume of upper basin box
    V2 = 77e9 - V1 # volume of deeper basin box
    # landward sill
    Hp = 49 # depth
    Bp = 1.4e3 # width
    Lp = 9e3 # length
    # landward basin
    V3 = 17e9 * 0.3 # volume of upper basin box
    V4 = 17e9 - V3 # volume of deeper basin box
    # constants
    beta = 7.7e-4
    g = 9.8
    Cd = 2.6e-3
    Socn = 32
    #
    dims = (H, B, L, Hp, Bp, Lp, V1, V2, V3, V4, beta, g, Cd, Socn)
    return dims

def fK(Cd, Ut, H):
    a0 = 0.028
    K = a0 * Cd * Ut * H
    Ks = K/2.2
    return K, Ks
    
def fSQin(H, B, K, Ks, L, Qr, Sout_l, Sin_s):
    c2 = g * beta * H * Sin_s
    A = H*B
    H2 = H**2
    L2 = L**2
    a = A*c2*H2/(4*48*Sin_s*K)
    b = 3*H2*a/(20*Ks*A)
    # find dsdx
    dsdx = (-L + np.sqrt(L2 + 8*b*(Sin_s-Sout_l)))/(4*b)
    dq = a * dsdx
    Qin = dq
    Qout = -(Qin + Qr)
    # adjust salinitied to conserve salt flux
    ds = b * dsdx**2
    eps = Qr*(Sin_s-Sout_l-2*ds)/(2*Qin+Qr)
    Sin_l = Sout_l + 2*ds - eps
    Sout_s = Sin_s - 2*ds - eps
    return Sin_l, Sout_s, Qin, Qout, dsdx
    
def advance_s(Qr, Ut, dims, S1, S2, S3, S4, dt, do_check=False, landward_basin=True):
    H, B, L, Hp, Bp, Lp, V1, V2, V3, V4, beta, g, Cd, Socn = dims
    # landward basin
    K, Ks = fK(Cd, Ut, Hp)
    Sout_l = S3
    Sin_s = S2
    Sin_lp, Sout_sp, Qinp, Qoutp, dsdxp = fSQin(Hp, Bp, K, Ks, Lp, Qr, Sout_l, Sin_s)
    S3n = S3 + dt*( S3*Qoutp/V1 + S4*Qinp/V1)
    S4n = S4 + dt*( Sin_lp*Qinp/V2 - S4*Qinp/V2)
    if landward_basin == False:
        S3n = 0
        S4n = 0
        Sout_sp = 0
        Qoutp = Qr
        Qinp = 0
    # seaward basin
    K, Ks = fK(Cd, Ut, H)
    Sout_l = S1
    Sin_s = Socn
    Sin_l, Sout_s, Qin, Qout, dsdx = fSQin(H, B, K, Ks, L, Qr, Sout_l, Sin_s)
    Q21 = Qin-Qinp
    if Q21 >= 0:
        F21 = Q21*S2
    else:
        F21 = Q21*S1
    S1n = S1 + dt*( S1*Qout/V1
        - Sout_sp*Qoutp/V1
        + F21/V1 )
    S2n = S2 + dt*( Sin_l*Qin/V2 - S2*Qin/V2)
    if do_check:
        # Checking on salt flux conservation at both ends
        # of the seaward strait
        # RESULT: it works perfectly now.
        Fl = Qin*Sin_l + Qout*Sout_l
        Fs = Qin*Sin_s + Qout*Sout_s
        print('F landward end = ' + str(Fl))
        print('F seaward end = ' + str(Fs))
    return S1n, S2n, S3n, S4n, Qin, Qinp

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

if __name__ == '__main__':
    # example use of the functions
    x = cubic_solver(1,2,3,4)
    print('The real root of the cubic is ' + str(the_root))

