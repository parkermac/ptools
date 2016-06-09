# -*- coding: utf-8 -*-
"""
Code for general plotting of functions.

Created on Sun Jan 24 07:34:10 2016

@author: PM5
"""

#%% IMPORTS

import matplotlib.pyplot as plt
import numpy as np

#%% CREATE DATA

# First make the x-values
x = np.linspace(-100, 100, 1000)

# Then define y = f(x)
# Example: y = 3*x**2 + 7*(x-1) + 5

y = (x-3) + (7/x) - 10

#%% PLOTTING

plt.close() # close old plots

plt.plot(x, y, '-r') # plot a red line

plt.grid() # add gridlines
plt.show() # make the plot visible

