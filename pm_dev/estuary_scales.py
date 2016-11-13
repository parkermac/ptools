#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 07:31:17 2016

@author: PM5

Calculate scales of an estuary.
"""

import numpy as np

g = 9.8
beta = 7.7e-4
omega = 1.4e-4
Cd = 3e-3
socn = 35

def get_M(g, beta, omega, Cd, Ut, socn, H):
    N0 = np.sqrt(beta*g*socn/H)
    M = np.sqrt( Cd*Ut**2/(omega*N0*H**2) )
    return M
    
def get_ubar(Qr, A):
    ubar = Qr/A
    return ubar
    
def get_c(g, beta, socn, H):
    c = np.sqrt(g*beta*socn*H)
    return c
    
def get_F(ubar, c):
    F = ubar/c
    return F

def get_FM(g, beta, omega, Cd, Ut, socn, H, Qr, A):   
    M = get_M(g, beta, omega, Cd, Ut, socn, H)
    ubar = get_ubar(Qr, A)
    c = get_c(g, beta, socn, H)
    F = get_F(ubar, c)
    print(('\n' + 5*'=' + ' %s ' + 5*'=') % (name))
    print('F = %8.1e, M = %8.1e' % (F, M))
    
 
H = 10
A = H * 10e3 / 2
Qr = 1500
    
name = 'aestus1, middle, spring'
Ut = 1.2
get_FM(g, beta, omega, Cd, Ut, socn, H, Qr, A)

name = 'aestus1, middle, neap'
Ut = 0.6
get_FM(g, beta, omega, Cd, Ut, socn, H, Qr, A)




