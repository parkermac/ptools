"""
Simple plotting program, that introduces basic matplotlib concepts.
"""

import matplotlib.pyplot as plt

plt.close()

fig1 = plt.figure()

ax1 = fig1.add_subplot(111)

x = 0
y = 0
ax1.text(x,y,'X',
fontsize=40,
horizontalalignment='center',
verticalalignment='center'
)

ax1.plot([-3, 3],[-1, -1],'-b',linewidth=3)
ax1.plot([-3, 3],[1, 1],'-b',linewidth=3)
ax1.plot([-1, -1],[-3, 3],'-b',linewidth=3)
ax1.plot([1, 1],[-3, 3],'-b',linewidth=3)

ax1.set_xlim((-5, 5))

ax1.set_ylim((-5, 5))

#ax1.grid()

plt.show()