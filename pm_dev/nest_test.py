#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 16:12:24 2016

@author: PM5

Code to test an offline two-way nesting scheme for LiveOcean.
"""

import numpy as np
import matplotlib.pyplot as plt

f = [0]
t = [0]

dt = .1
T = 10
om = 2*np.pi/T

while t[-1] < 100:
    if t[-1] == 0:
        f.append(f[-1] + dt*np.cos(om*t[-1]))
    else:
        f.append(f[-1] + dt*(1*(f[-1]-f[-2]) + np.cos(om*t[-1])))
    t.append(t[-1] + dt)

f = np.array(f)
t = np.array(t)

fa = (1/om)*np.sin(om*t)
  
plt.close()

plt.plot(t,f,'-r', t,fa,'-b')

plt.show()

