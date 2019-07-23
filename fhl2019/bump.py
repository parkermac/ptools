"""
Code to explore shallow water hydraulics.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

# set parameters
g = 9.8 # gravity [m/s2]
H0 = 10 # far field channel depth [m]
U_dict = {0:3, 1:20} # Far field velocities [m/s] [Subcritical, Supercritical]
L = 1 # channel length [arbitrary units]
Lb = L/10 # bump length
x = np.linspace(-L/2, L/2, 500) # x-coordinate [m]
NX = len(x)

def cubic_solver(a,b,c,d):
    """
    This gives the real roots of a cubic polynomial of the form
    a*x**3 + b*x**2 + c*x + d = 0.
    """
    b3a = b/(3*a)
    b3a2 = b3a * b3a
    b3a3 = b3a * b3a * b3a
    p = 3*b3a2 - 2*(b/a)*b3a + c/a
    q = -b3a3 + (b/a)*b3a2 - (c/a)*b3a + d/a
    D = (p/3)**3 + (q/2)**2
    signD = np.sign(D)
    if D > 0:
        # The equation has one real root.
        u_inn = -q/2 + np.sqrt(D)
        u = np.sign(u_inn)*np.abs(np.complex(u_inn)**(1/3))
        v_inn = -q/2 - np.sqrt(D)
        v = np.sign(v_inn)*np.abs(np.complex(v_inn)**(1/3))
    else:
        # The equation has three distinct real roots,
        # unless D=0 in which case at least two are equal.
        u_inn = -q/2 + np.sqrt(np.complex(D))
        u = np.complex(u_inn)**(1/3)
        v_inn = -q/2 - np.sqrt(np.complex(D))
        v = np.complex(v_inn)**(1/3)
    # there are three roots (note 1j is python's ways of saying i)
    F1 = u + v 
    F2 = -(u + v)/2 + (u - v)*1j*np.sqrt(3)/2
    F3 = -(u + v)/2 - (u - v)*1j*np.sqrt(3)/2
    rr = np.nan * np.ones(3)
    rr[0] = np.real(F1 - b3a)
    if D <=0:
        rr[1] = np.real(F2 - b3a)
        rr[2] = np.real(F3 - b3a)
    
    return rr, signD

plt.close('all')
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(15,8))

counter = 0
for delta_b in np.linspace(.3, 5, 6): # bump height [m]

    for iu in range(2):
        U = U_dict[iu]
        
        Q = U * H0 # transport per unit width
        # bathymetry (defined positive up, and relative to z=-H0)
        delta = delta_b * np.exp(-(x*x)/(Lb*Lb))
        eta = np.nan * np.ones_like(x)
        for ii in range(NX):
            # coefficients of the cubic
            AA = g
            BB = -(U**2)/2 - g*H0 + g*delta[ii]
            CC = 0
            DD = (Q**2)/2
            rr, signD = cubic_solver(AA, BB, CC, DD)
            # figure out which root is the right one
            if signD > 0: # one real root
                ir = 0
            else: # three real roots
                ir = (np.abs(rr-H0)).argmin()
            hh = rr[ir]
            if hh > 0:
                eta[ii] = hh - H0 + delta[ii]
            else:
                eta[ii] = np.nan

        h = H0 + eta - delta
        cc = np.sqrt(g*h)
        u = Q/h
        F = u/cc

        # plotting
        alpha = .1 + .1 * counter
        
        if counter == 0:
            #
            ax = axes[0,iu]
            ax.plot(x, -H0+delta, '-k', linewidth=2, alpha=alpha)
            ax.plot(x, eta, '-c', linewidth=3, alpha=alpha)
            ax.text(.05, .7, 'Free surface [m]', color='c',
                transform=ax.transAxes, fontweight='bold')
            ax.text(.05, .3, 'Bottom [m]', color='k',
                transform=ax.transAxes, fontweight='bold')
            ax.set_xlim(-L/2, L/2)
            if iu == 0:
                ax.set_title('Subcritical')
            elif iu == 1:
                ax.set_title('Supercritical')
            #
            ax = axes[1,iu]
            ax.plot(x, F, '-r', linewidth=3, alpha=alpha)
            ax.text(.05, .5, 'Froude Number', color='r',
                transform=ax.transAxes, fontweight='bold')
            ax.set_xlim(-L/2, L/2)
            #
            ax = axes[2,iu]
            ax.plot(x, u, '-b', linewidth=3, alpha=alpha)
            ax.text(.05, .5, 'U [m/s]', color='b',
                transform=ax.transAxes, fontweight='bold')
            ax.set_xlim(-L/2, L/2)
            ax.set_xlabel('Along Channel Distance')
        else:
            ax = axes[0,iu]
            ax.plot(x, -H0+delta, '-k', linewidth=2, alpha=alpha)
            ax.plot(x, eta, '-c', linewidth=3, alpha=alpha)
            #
            ax = axes[1,iu]
            ax.plot(x, F, '-r', linewidth=3, alpha=alpha)
            #
            ax = axes[2,iu]
            ax.plot(x, u, '-b', linewidth=3, alpha=alpha)
        
    counter += 1

plt.show()

