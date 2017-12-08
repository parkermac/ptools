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