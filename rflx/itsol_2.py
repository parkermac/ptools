"""
Code to iteratively solve for the flux fraction coefficients
at a channel segment.

Inovlves smarter initial guesses.
"""

import numpy as np

# Specify the inflow (i) and outflow(o) transport and salinity.
# The index into each array corresponds to which section is it.

case = 2

if case==2:
    # For this example we expect:
    #   Abest[0,0] = alpha34 = .4
    #   Abest[1,1] = alpha21 = .6
    Qi = np.array([25, 25])
    Qo = np.array([20, 30])
    Si = np.array([20,30])
    So = np.array([25, 25])
elif case==3:
    # Triple junction
    Qi = np.array([25, 25, 2])
    Qo = np.array([20, 30, 2])
    Si = np.array([20,30, 27])
    So = np.array([25, 25, 22])
elif case==4:
    # Quadruple junction
    Qi = np.array([25, 25, 2, 5])
    Qo = np.array([20, 30, 2, 5])
    Si = np.array([20,30, 27, 24])
    So = np.array([25, 25, 22, 26])

# salt fluxes
Fi = Qi * Si
Fo = Qo * So

def get_a0(Qo, Fo):
    # Create initial guess for transport fraction vector.
    ns = len(Qo) # number of sections
    aq = np.nan * np.ones((ns*ns))
    af = np.nan * np.ones((ns*ns))
    for oo in range(ns):
        for ii in range(ns):
            aq[oo*ns + ii] = Qo[ii] / Qo.sum()
            af[oo*ns + ii] = Fo[ii] / Fo.sum()
    a = (aq + af)/2
    return a
    
def get_err(a, qi, qo, fi, fo):
    # calculate the error
    Qe = Qi.dot(a) - Qo
    Fe = Fi.dot(a) - Fo
    return Qe, Fe
    
a = get_a0(Qo, Fo)
Qe, Fe = get_err(a, Qi, Qo, Fi, Fo)
    
print('\na:')
print(a)
print('\nqe:')
print(Qe)
print('\nfe:')
print(Fe)
