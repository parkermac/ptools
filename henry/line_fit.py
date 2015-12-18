"""
Simple plotting code, with a best fit line
"""

# define the data
import numpy as np
x = np.array([14501, 127, 248, 1763, 447])
y = np.array([0, 4601, 4253, 1582, 3734])

# Fit the data with a line
#
# polyfit will fit any order of polynomial, and by making the fit_order = 1
# we are telling it to fit a first order polynomial, which is
# a straight line.  If we had  used fit_order = 2 it would fit a parabola.
#
fit_order = 1
p = np.polyfit(x, y, fit_order)
# The returned array p is the coefficients of the line,
# arranged so that the first item p[0] is the coefficient of the highest-order
# term.  For a linear fit the highest order term is the one with "x" in it, so
# p[0] would be the slope, and p[1] would be the y-intercept.

# create a line for the fit:
#
# first make the x axis for the best fit line
# (call it xx, so it is not confused with the data)
xx = np.linspace(x.min(), x.max(), 100)
# then the y axis
yy = np.polyval(p, xx)
#yy = p[0]*xx + p[1] # like y = mx * b

# plotting
import matplotlib.pyplot as plt
plt.close() # get rid of old plots
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)

ax.plot(x, y, '*r', markersize=20) # plot the raw data
ax.plot(xx, yy, '-g', linewidth=3) # plot the fit line
ax.grid() # add gridlines for fun

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('Raw Data and Best Fit Line')

# put the fit equation on the plot
ax.text(0.9, 0.9,
    '$y=%0.2fx+%0.2f$' % (p[0], p[1]), # the $ makes it use LaTeX formatting
    horizontalalignment='right',
    fontsize=20,
    color='g',
    transform=ax.transAxes)
    
plt.show() # make sure the plot appears