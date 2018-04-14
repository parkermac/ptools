#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 16:56:31 2018

@author: pm7

Code to test the streamplot function.

"""

import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-10, 10, num=100)
y = x

X, Y = np.meshgrid(x, y)


u = X**2 + np.sin(Y)
v = -X + Y/2

z = np.sqrt(u**2 + v**2)

plt.close('all')
fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(111)

ax.streamplot(x, y, u,v, density=2,
              color=z, linewidth=z/25,
              arrowstyle='-', cmap='rainbow')

plt.show()