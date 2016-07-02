# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 14:09:25 2016

@author: PM5

Code to test speed of numpy.
"""

import numpy as np

N = 1000
M = 1000

a = np.random.randn(N, N)

aa = np.sin(a)

import time
tt0 = time.time()

for ii in range(M):
    aa = np.sin(aa)

print('Took %0.1f seconds' % (time.time() - tt0))