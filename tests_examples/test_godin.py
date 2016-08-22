# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 11:54:40 2016

@author: PM5

Test of the rewritten Hanning and Godin filtering functions.

RESULT: they appear to work correctly, and now the Hanning filter can accept
even or odd lengths.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
alp = os.path.abspath('/Users/PM5/Documents/LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
from importlib import reload
import zfun
reload(zfun)

# make data
t = np.arange(24*100)
data = 1*np.cos(2*np.pi*t/60) + 1*np.cos(2*np.pi*t/12.42) + .5*np.cos(2*np.pi*t/24)

# godin filtering by hand
filt = zfun.godin_shape()
filt = filt / filt.sum()
n = np.ceil(len(filt)/2).astype(int)
smooth = np.convolve(data, filt, mode = 'same')
smooth[:n] = np.nan
smooth[-n:] = np.nan

# and test out the zfun versions
smooth_h40 = zfun.filt_hanning(data)
smooth_h71 = zfun.filt_hanning(data, n=71)
smooth_g = zfun.filt_godin(data)

plt.close()
plt.plot(t, data, '-k', t, smooth, '-r')
plt.plot(t, smooth_h40, '-g', t, smooth_h71, '-b', t, smooth_g, '-m')

plt.show()