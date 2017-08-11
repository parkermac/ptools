"""
Code to test setting default font sizes and etc. in 
matplotlib figures

"""

import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 10, 100)
y = np.sin(x)

# PLOTTING

plt.close('all')

# SET RC
fs1 = 16
fs2 = 20
lw1 = 3
lw2 = 5
plt.rc('xtick', labelsize=fs1)
plt.rc('ytick', labelsize=fs1)
plt.rc('xtick.major', size=10, pad=5, width=lw1)
plt.rc('ytick.major', size=10, pad=5, width=lw1)
plt.rc('axes', lw=lw1)
plt.rc('lines', lw=lw2)
plt.rc('font', size=fs2)
plt.rc('grid', color='g', ls='-', lw=lw1, alpha=.3)

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
ax.plot(x, y, '-b')
ax.grid()
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_xlim(0, 10)
ax.set_ylim(-1.5, 1.5)

# RC CLEANUP
plt.tight_layout()
plt.rcdefaults()

plt.show()