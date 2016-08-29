# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 09:58:03 2016

@author: PM5

Code to test the statistical functions
"""

import numpy as np
import matplotlib.pyplot as plt

from importlib import reload
import statfun as sfun
reload(sfun)

N = 500
x = np.linspace(0, 100, N)

a = 1
b = 0.5
c = 0.5

def make_y(a, b, c, N):
    y = a * np.cos(x/5) + b * np.cos(x/8) + c * np.random.randn(N)
    return y

y1 = make_y(a, b, c, N)
y2 = make_y(a, b, c, N)

#%% plotting

plt.close()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x, y1, '-r', x, y2, '-b')

cc = sfun.cross_correlation(y1, y2)

ax.set_title('cc = %0.2f' % (cc))

plt.show()

