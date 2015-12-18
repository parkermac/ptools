"""
Plotting Series
"""

import numpy as np
import matplotlib.pyplot as plt

# set up the constants that define the sequence or series
a1 = .5
r = .9

# this formula gives the exact infinite sum of the series
Sinf = a1 / (1-r)

# this calculates the sum for each step along the way

n = 100 # max number of steps
ivec = range(n)
iivec = range(1,n+1)

s = np.zeros(n) # terms in the sequence (not summed yet)

s[0] = a1 # set the first value to be a1

# now make each of the next terms in the sequence
for ii in ivec[1:]:
    s[ii] = s[ii-1] * r
    #print('ii = ' + str(ii) + ' s = ' + str(s[ii]))

# now add them up to make the series
SS = np.cumsum(s)

plt.close()
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

ax.plot(iivec, SS, '-*r', linewidth=2)

SSinf = Sinf * np.ones(n)

ax.plot(iivec, SSinf, '-b', linewidth=2)

plt.show()